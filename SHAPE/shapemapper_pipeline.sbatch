#!/bin/bash

# ===================================================================================
# Description:
#   This script orchestrates a full SHAPE-Mapper analysis pipeline, from running
#   SHAPE-Mapper to generating final structure models with RNAfold. It is designed
#   to be submitted to a Slurm cluster as a job array, with each job processing
#   one full sample (e.g., HEK293T, which includes Rep1 and Rep2).
#
#   The pipeline reads a sample map file to get sample names and FASTQ paths.
#
# Arguments:
#   --sample-map <path>         Path to the TSV file mapping sample names to FASTQ files. (Required)
#   --output-dir <path>         Base directory for all pipeline outputs. (Required)
#   --ref-fasta-dir <path>      Directory containing FASTA chunks for SHAPE-Mapper. (Required)
#   --ref-rnafold-fasta <path>  Single FASTA file for RNAfold. (Required)
#   --num-groups <int>          Total number of unique sample groups to process. (Required)
#   --group-col <int>           Column number for the group name in the sample map (Default: 1)
#   -h, --help                  Display this help message and exit.
#
# Sample Map Format (TSV):
#   A header is required. The default columns are:
#   1: GroupName (shared ID for all replicates of a sample)
#   2: SampleName (unique ID for each replicate)
#   3: UntreatedFastq (path)
#   4: ModifiedFastq (path)
#
# SLURM Usage:
#   sbatch --array=1-N --mail-user=rodell@stanford.edu master_pipeline.sh \
#       --sample-map /path/to/sample_map.tsv \
#       --output-dir /path/to/pipeline_output \
#       --ref-fasta-dir /path/to/fasta_chunks \
#       --ref-rnafold-fasta /path/to/rnafold.fasta \
#       --num-groups N
#   (Where N is the number of unique groups, e.g., 2 for HEK293T and HepG2)
#
# ===================================================================================

# --- SLURM Settings ---
#SBATCH --job-name=shapemapper_pipeline_%j
#SBATCH --output=shapemapper_pipeline_%j.out
#SBATCH --error=shapemapper_pipeline_%j.err
#SBATCH --time=4:00:00
#SBATCH -p normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=16GB
#SBATCH --mail-type=BEGIN,END,FAIL

set -e # Exit immediately if a command exits with a non-zero status.
set -o pipefail # Fail a pipeline if any command in it fails.

# --- Define a function for displaying usage information ---
usage() {
    echo "Usage: $0 --sample-map <path> --output-dir <path> --ref-fasta-dir <path> --ref-rnafold-fasta <path> --num-samples <int>"
    echo ""
    echo "Required Arguments:"
    echo "  --sample-map <path>         Path to the sample map TSV file."
    echo "  --output-dir <path>         Base directory for all pipeline outputs."
    echo "  --ref-fasta-dir <path>      Directory with FASTA chunks for SHAPE-Mapper."
    echo "  --ref-rnafold-fasta <path>  FASTA file for RNAfold."
    echo "  --num-groups <int>          Total number of unique sample groups for the array."
    echo ""
    echo "Optional Arguments:"
    echo "  --group-col <int>             Column number for the GroupName in the sample map (Default: 1)."
    echo "  -h, --help                    Display this help message and exit."
    exit 1
}

# --- 1. ARGUMENT PARSING AND VALIDATION ---
SAMPLE_MAP_FILE=""
BASE_PIPELINE_DIR=""
REF_FASTA_DIR=""
REF_RNAFOLD_FASTA=""
NUM_GROUPS=""
GROUP_COL=1 # Default column for GroupName

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --sample-map) SAMPLE_MAP_FILE="$2"; shift ;;
        --output-dir) BASE_PIPELINE_DIR="$2"; shift ;;
        --ref-fasta-dir) REF_FASTA_DIR="$2"; shift ;;
        --ref-rnafold-fasta) REF_RNAFOLD_FASTA="$2"; shift ;;
        --num-groups) NUM_GROUPS="$2"; shift ;;
        --group-col) GROUP_COL="$2"; shift ;;
        -h|--help) usage ;;
        *) echo "ERROR: Unknown parameter passed: $1"; usage ;;
    esac
    shift
done

# Check for missing required arguments
if [[ -z "$SAMPLE_MAP_FILE" || -z "$BASE_PIPELINE_DIR" || -z "$REF_FASTA_DIR" || -z "$REF_RNAFOLD_FASTA" || -z "$NUM_GROUPS" ]]; then
    echo "ERROR: One or more required arguments are missing."
    usage
fi

# Validate paths and Slurm environment
if [[ ! -f "$SAMPLE_MAP_FILE" ]]; then echo "ERROR: Sample map file not found: ${SAMPLE_MAP_FILE}"; exit 1; fi
if [[ ! -d "$REF_FASTA_DIR" ]]; then echo "ERROR: Reference FASTA directory not found: ${REF_FASTA_DIR}"; exit 1; fi
if [[ ! -f "$REF_RNAFOLD_FASTA" ]]; then echo "ERROR: RNAfold reference FASTA not found: ${REF_RNAFOLD_FASTA}"; exit 1; fi
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then echo "ERROR: This script must be run as a Slurm job array."; exit 1; fi
if [ "$SLURM_ARRAY_TASK_COUNT" -ne "$NUM_GROUPS" ]; then echo "Warning: --num-groups (${NUM_GROUPS}) does not match Slurm array count (${SLURM_ARRAY_TASK_COUNT})."; fi

# --- 2. DEFINE PATHS AND VARIABLES ---
# --- Script Directories ---
RNASTRUCTURE_DIR="$HOME/RNAstructure"
SHAPE_SCRIPT_DIR="${RNASTRUCTURE_DIR}/SHAPE"

# --- Dynamic Paths based on Sample ---
UNIQUE_GROUPS=($(tail -n +2 "$SAMPLE_MAP_FILE" | cut -f${GROUP_COL} | sort -u))
GROUP_NAME=${UNIQUE_GROUPS[$((SLURM_ARRAY_TASK_ID - 1))]}

if [ -z "$GROUP_NAME" ]; then
    echo "ERROR: Array task ID ${SLURM_ARRAY_TASK_ID} is out of bounds for the sample list."
    exit 1
fi


# --- Output Directories for this specific sample ---
GROUP_OUT_DIR="${BASE_PIPELINE_DIR}/${GROUP_NAME}"
SHAPEMAPPER_OUT_DIR="${GROUP_OUT_DIR}/shapemapper"
RNAFOLD_OUT_DIR="${GROUP_OUT_DIR}/rnafold"
PLOTS_DIR="${GROUP_OUT_DIR}/plots"
SHAPE_DIR="${GROUP_OUT_DIR}/shape_files"
LOG_DIR="${GROUP_OUT_DIR}/logs"

mkdir -p "$SHAPEMAPPER_OUT_DIR" "$RNAFOLD_OUT_DIR" "$PLOTS_DIR" "$SHAPE_DIR" "$LOG_DIR"

# --- Setup Logging ---
LOG_FILE="${LOG_DIR}/${GROUP_NAME}_pipeline.log"
exec &> >(tee -a "${LOG_FILE}") # Redirect all output to log and console

echo "================================================================="
echo "STARTING PIPELINE for SAMPLE: ${GROUP_NAME}"
echo "Array Job ID: ${SLURM_ARRAY_TASK_ID}"
echo "Log file: ${LOG_FILE}"
echo "================================================================="

# --- 3. FUNCTION TO PROCESS A SINGLE REPLICATE (Steps 1-4a) ---
process_replicate() {
    local rep_name="$1"
    echo "--- [START] Processing Replicate: ${rep_name} ---"

    map_line=$(awk -F'\t' -v name="$rep_name" 'NR > 1 && $2 == name' "$SAMPLE_MAP_FILE")
    untreated_fq=$(echo "$map_line" | cut -f3)
    modified_fq=$(echo "$map_line" | cut -f4)

    if [[ -z "$map_line" || ! -f "$untreated_fq" || ! -f "$modified_fq" ]]; then
        echo "ERROR: FASTQ files for '${rep_name}' not found or line could not be parsed from sample map."
        echo "DEBUG: Extracted Untreated Path: '${untreated_fq}'"
        echo "DEBUG: Extracted Modified Path: '${modified_fq}'"
        exit 1
    fi
    echo "Found FASTQs: Untreated=${untreated_fq}, Modified=${modified_fq}"

    echo "[${rep_name}] Step 1: Running SHAPE-Mapper..."
    bash "${SHAPE_SCRIPT_DIR}/run_shapemapper.sh" "$modified_fq" "$untreated_fq" "$rep_name" "$SHAPEMAPPER_OUT_DIR" "$REF_FASTA_DIR"

    echo "[${rep_name}] Step 2: Combining SHAPE profiles..."
    Rscript "${SHAPE_SCRIPT_DIR}/02_SHAPEprofile_combine.R" --dir "$SHAPEMAPPER_OUT_DIR" --pattern "${rep_name}_chunk_" --outfile "${GROUP_OUT_DIR}/${rep_name}_profile.txt"

    echo "[${rep_name}] Step 3: Extracting poor quality RNA list..."
    # bash "${SHAPE_SCRIPT_DIR}/03_extract_poor_quality_rnas.sh" --log-file "${SHAPEMAPPER_OUT_DIR}/logs/${rep_name}_batch_run.log" --outfile "${GROUP_OUT_DIR}/${rep_name}_poor_quality_rna.txt"
    grep "possible poor quality reactivity profiles" "${SHAPEMAPPER_OUT_DIR}/logs/${rep_name}_chunk_"*.log | \
    awk -F':' '{sub(/^[ ]+/, "", $NF); print $NF}' | \
    awk '{gsub(/, /,"\n"); print}' | \
    sort -u > "${GROUP_OUT_DIR}/${rep_name}_poor_quality_rna.txt"
    
    echo "[${rep_name}] Step 4: Annotating profiles..."
    Rscript "${SHAPE_SCRIPT_DIR}/04_SHAPE_annotate_profiles.R" --profile "${GROUP_OUT_DIR}/${rep_name}_profile.txt" --list "${GROUP_OUT_DIR}/${rep_name}_poor_quality_rna.txt" --outfile "${GROUP_OUT_DIR}/${rep_name}_profile_annotated.txt"

    echo "--- [DONE] Processing Replicate: ${rep_name} ---"
}

# --- 4. MAIN EXECUTION BLOCK ---
# Load modules once for the entire job
ml biology viennarna R

# Find all replicate sample names for the current group
REPLICATE_NAMES=($(awk -F'\t' -v group="$GROUP_NAME" -v group_col="$GROUP_COL" 'NR > 1 && $group_col == group {print $2}' "$SAMPLE_MAP_FILE"))
NUM_REPS=${#REPLICATE_NAMES[@]}

echo "Found ${NUM_REPS} sample(s) for group '${GROUP_NAME}': ${REPLICATE_NAMES[*]}"

if [[ "$NUM_REPS" -eq 0 ]]; then
    echo "ERROR: No samples found for group '${GROUP_NAME}'. Check your sample map."
    exit 1
fi

if [[ "$NUM_REPS" -eq 2 ]]; then
    # --- Replicated Workflow (2 Replicates) ---
    echo "Detected 2 replicates. Running replicated workflow for ${GROUP_NAME}."
    
    REP1_NAME=${REPLICATE_NAMES[0]}
    REP2_NAME=${REPLICATE_NAMES[1]}

    process_replicate "$REP1_NAME" &
    PID_REP1=$!
    process_replicate "$REP2_NAME" &
    PID_REP2=$!
    
    echo "Waiting for replicate processing (PIDs: ${PID_REP1}, ${PID_REP2})..."
    wait $PID_REP1 $PID_REP2
    echo "All replicate processing finished."

    echo "Step 4b: Comparing replicate reactivities..."
    Rscript "${SHAPE_SCRIPT_DIR}/04b_compare_replicates.R" --rep1 "${GROUP_OUT_DIR}/${REP1_NAME}_profile_annotated.txt" --rep2 "${GROUP_OUT_DIR}/${REP2_NAME}_profile_annotated.txt" --out-prefix "${PLOTS_DIR}/${GROUP_NAME}"
    
    echo "Step 5: Averaging reactivities across replicates..."
    Rscript "${SHAPE_SCRIPT_DIR}/05_average_reactivities.R" --infiles "${GROUP_OUT_DIR}/${REP1_NAME}_profile_annotated.txt" "${GROUP_OUT_DIR}/${REP2_NAME}_profile_annotated.txt" --outdir "$SHAPE_DIR"

elif [[ "$NUM_REPS" -eq 1 ]]; then
    # --- Single Sample Workflow ---
    echo "Detected 1 sample. Running single-sample workflow for ${GROUP_NAME}."
    
    REP1_NAME=${REPLICATE_NAMES[0]}
    process_replicate "$REP1_NAME"
    
    echo "Skipping Step 4b (Replicate Comparison)."
    echo "Step 5 (single): Generating .shape files from single profile..."
    Rscript "${SHAPE_SCRIPT_DIR}/05_average_reactivities.R" --infiles "${GROUP_OUT_DIR}/${REP1_NAME}_profile_annotated.txt" --outdir "$SHAPE_DIR"

else
    echo "ERROR: This script currently supports 1 or 2 replicates per group, but found ${NUM_REPS} for '${GROUP_NAME}'."
    echo "Please adjust the script or your sample map."
    exit 1
fi


# --- Step 6: SHAPE-informed RNAfold ---
echo "Step 6: Running SHAPE-informed RNAfold for ${GROUP_NAME}..."
bash "${SHAPE_SCRIPT_DIR}/06_SHAPEinformed_RNAfold.sh" -f "$REF_RNAFOLD_FASTA" -s "$SHAPE_DIR" -o "$RNAFOLD_OUT_DIR"

# --- Step 7: Extract RNAfold Summaries ---
echo "Step 7a: Extracting pairing probabilities..."
Rscript "${RNASTRUCTURE_DIR}/extract_pairing_prob.R" --indir "$RNAFOLD_OUT_DIR" --outfile "${GROUP_OUT_DIR}/${GROUP_NAME}_pairing_probabilities.csv" --max-len 130 # Adjust if needed

echo "Step 7b: Extracting fold summaries..."
Rscript "${RNASTRUCTURE_DIR}/fold_summary.R" --indir "$RNAFOLD_OUT_DIR" --outfile "${GROUP_OUT_DIR}/${GROUP_NAME}_rnafold_summary.csv"

echo "================================================================="
echo "PIPELINE COMPLETE for SAMPLE: ${GROUP_NAME}"
echo "================================================================="