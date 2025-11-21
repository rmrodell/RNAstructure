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
#   --num-samples <int>         Total number of unique samples to process (for array validation). (Required)
#   -h, --help                  Display this help message and exit.
#
# SLURM Usage:
#   sbatch --array=1-N --mail-user=rodell@stanford.edu master_pipeline.sh \
#       --sample-map /path/to/sample_map.tsv \
#       --output-dir /path/to/pipeline_output \
#       --ref-fasta-dir /path/to/fasta_chunks \
#       --ref-rnafold-fasta /path/to/rnafold.fasta \
#       --num-samples N
#   (Where N is the number of unique samples, e.g., 2 for HEK293T and HepG2)
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
    echo "  --num-samples <int>         Total number of unique samples for the array."
    echo ""
    echo "Optional Arguments:"
    echo "  -h, --help                    Display this help message and exit."
    exit 1
}

# --- 1. ARGUMENT PARSING AND VALIDATION ---
SAMPLE_MAP_FILE=""
BASE_PIPELINE_DIR=""
REF_FASTA_DIR=""
REF_RNAFOLD_FASTA=""
NUM_SAMPLES=""

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --sample-map) SAMPLE_MAP_FILE="$2"; shift ;;
        --output-dir) BASE_PIPELINE_DIR="$2"; shift ;;
        --ref-fasta-dir) REF_FASTA_DIR="$2"; shift ;;
        --ref-rnafold-fasta) REF_RNAFOLD_FASTA="$2"; shift ;;
        --num-samples) NUM_SAMPLES="$2"; shift ;;
        -h|--help) usage ;;
        *) echo "ERROR: Unknown parameter passed: $1"; usage ;;
    esac
    shift
done

# Check for missing required arguments
if [[ -z "$SAMPLE_MAP_FILE" || -z "$BASE_PIPELINE_DIR" || -z "$REF_FASTA_DIR" || -z "$REF_RNAFOLD_FASTA" || -z "$NUM_SAMPLES" ]]; then
    echo "ERROR: One or more required arguments are missing."
    usage
fi

# Validate paths and Slurm environment
if [[ ! -f "$SAMPLE_MAP_FILE" ]]; then echo "ERROR: Sample map file not found: ${SAMPLE_MAP_FILE}"; exit 1; fi
if [[ ! -d "$REF_FASTA_DIR" ]]; then echo "ERROR: Reference FASTA directory not found: ${REF_FASTA_DIR}"; exit 1; fi
if [[ ! -f "$REF_RNAFOLD_FASTA" ]]; then echo "ERROR: RNAfold reference FASTA not found: ${REF_RNAFOLD_FASTA}"; exit 1; fi
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then echo "ERROR: This script must be run as a Slurm job array."; exit 1; fi
if [ "$SLURM_ARRAY_TASK_COUNT" -ne "$NUM_SAMPLES" ]; then echo "Warning: --num-samples (${NUM_SAMPLES}) does not match Slurm array count (${SLURM_ARRAY_TASK_COUNT})."; fi

# --- 2. DEFINE PATHS AND VARIABLES ---
# --- Script Directories ---
RNASTRUCTURE_DIR="$HOME/RNAstructure"
SHAPE_SCRIPT_DIR="${RNASTRUCTURE_DIR}/SHAPE"

# --- Dynamic Paths based on Sample ---
UNIQUE_SAMPLES=($(tail -n +2 "$SAMPLE_MAP_FILE" | cut -f1 | sed 's/_Rep[12]$//' | sort -u))
SAMPLE_NAME=${UNIQUE_SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}

if [ -z "$SAMPLE_NAME" ]; then
    echo "ERROR: Array task ID ${SLURM_ARRAY_TASK_ID} is out of bounds for the sample list."
    exit 1
fi

SAMPLE_REP1="${SAMPLE_NAME}_Rep1"
SAMPLE_REP2="${SAMPLE_NAME}_Rep2"

# --- Output Directories for this specific sample ---
SAMPLE_OUT_DIR="${BASE_PIPELINE_DIR}/${SAMPLE_NAME}"
SHAPEMAPPER_OUT_DIR="${SAMPLE_OUT_DIR}/shapemapper"
RNAFOLD_OUT_DIR="${SAMPLE_OUT_DIR}/rnafold"
PLOTS_DIR="${SAMPLE_OUT_DIR}/plots"
SHAPE_DIR="${SAMPLE_OUT_DIR}/shape_files"
LOG_DIR="${SAMPLE_OUT_DIR}/logs"

mkdir -p "$SHAPEMAPPER_OUT_DIR" "$RNAFOLD_OUT_DIR" "$PLOTS_DIR" "$SHAPE_DIR" "$LOG_DIR"

# --- Setup Logging ---
LOG_FILE="${LOG_DIR}/${SAMPLE_NAME}_pipeline.log"
exec &> >(tee -a "${LOG_FILE}") # Redirect all output to log and console

echo "================================================================="
echo "STARTING PIPELINE for SAMPLE: ${SAMPLE_NAME}"
echo "Array Job ID: ${SLURM_ARRAY_TASK_ID}"
echo "Log file: ${LOG_FILE}"
echo "================================================================="

# --- 3. FUNCTION TO PROCESS A SINGLE REPLICATE (Steps 1-4a) ---
process_replicate() {
    local rep_name="$1"
    echo "--- [START] Processing Replicate: ${rep_name} ---"

    map_line=$(grep "^${rep_name}\b" "$SAMPLE_MAP_FILE")
    untreated_fq=$(echo "$map_line" | cut -f2)
    modified_fq=$(echo "$map_line" | cut -f3)

    if [[ ! -f "$untreated_fq" || ! -f "$modified_fq" ]]; then echo "ERROR: FASTQ files for ${rep_name} not found."; exit 1; fi
    echo "Found FASTQs: Untreated=${untreated_fq}, Modified=${modified_fq}"

    echo "[${rep_name}] Step 1: Running SHAPE-Mapper..."
    bash "${SHAPE_SCRIPT_DIR}/run_shapemapper.sh" "$modified_fq" "$untreated_fq" "$rep_name" "$SHAPEMAPPER_OUT_DIR"

    echo "[${rep_name}] Step 2: Combining SHAPE profiles..."
    Rscript "${SHAPE_SCRIPT_DIR}/02_SHAPEprofile_combine.R" --dir "$SHAPEMAPPER_OUT_DIR" --pattern "${rep_name}_chunk_" --outfile "${SAMPLE_OUT_DIR}/${rep_name}_profile.txt"

    echo "[${rep_name}] Step 3: Extracting poor quality RNA list..."
    # bash "${SHAPE_SCRIPT_DIR}/03_extract_poor_quality_rnas.sh" --log-file "${SHAPEMAPPER_OUT_DIR}/logs/${rep_name}_batch_run.log" --outfile "${SAMPLE_OUT_DIR}/${rep_name}_poor_quality_rna.txt"
    grep "possible poor quality reactivity profiles" "${SHAPEMAPPER_OUT_DIR}/logs/${rep_name}_chunk_"*.log | \
    awk -F':' '{sub(/^[ ]+/, "", $NF); print $NF}' | \
    awk '{gsub(/, /,"\n"); print}' | \
    sort -u > "${SAMPLE_OUT_DIR}/${rep_name}_poor_quality_rna.txt"
    
    echo "[${rep_name}] Step 4: Annotating profiles..."
    Rscript "${SHAPE_SCRIPT_DIR}/04_SHAPE_annotate_profiles.R" --profile "${SAMPLE_OUT_DIR}/${rep_name}_profile.txt" --list "${SAMPLE_OUT_DIR}/${rep_name}_poor_quality_rna.txt" --outfile "${SAMPLE_OUT_DIR}/${rep_name}_profile_annotated.txt"

    echo "--- [DONE] Processing Replicate: ${rep_name} ---"
}

# --- 4. MAIN EXECUTION BLOCK ---
# Load modules once for the entire job
ml biology viennarna R

# Detect if replicates exist
REP1_EXISTS=$(grep -c "^${SAMPLE_NAME}_Rep1\b" "$SAMPLE_MAP_FILE")
REP2_EXISTS=$(grep -c "^${SAMPLE_NAME}_Rep2\b" "$SAMPLE_MAP_FILE")

if [[ "$REP1_EXISTS" -eq 1 && "$REP2_EXISTS" -eq 1 ]]; then
    # --- Replicated Workflow ---
    echo "Detected Replicates: Running replicated workflow for ${SAMPLE_NAME}"
    
    SAMPLE_REP1="${SAMPLE_NAME}_Rep1"
    SAMPLE_REP2="${SAMPLE_NAME}_Rep2"

    process_replicate "$SAMPLE_REP1" &
    PID_REP1=$!
    process_replicate "$SAMPLE_REP2" &
    PID_REP2=$!
    
    echo "Waiting for replicate processing (PIDs: ${PID_REP1}, ${PID_REP2})..."
    wait $PID_REP1 $PID_REP2
    echo "All replicate processing finished."

    echo "Step 4b: Comparing replicate reactivities..."
    Rscript "${SHAPE_SCRIPT_DIR}/04b_compare_replicates.R" --rep1 "${SAMPLE_OUT_DIR}/${SAMPLE_REP1}_profile_annotated.txt" --rep2 "${SAMPLE_OUT_DIR}/${SAMPLE_REP2}_profile_annotated.txt" --out-prefix "${PLOTS_DIR}/${SAMPLE_NAME}"
    
    echo "Step 5: Averaging reactivities across replicates..."
    Rscript "${SHAPE_SCRIPT_DIR}/05_average_reactivities.R" --infiles "${SAMPLE_OUT_DIR}/${SAMPLE_REP1}_profile_annotated.txt" "${SAMPLE_OUT_DIR}/${SAMPLE_REP2}_profile_annotated.txt" --outdir "$SHAPE_DIR"

else
    # --- Single Sample Workflow ---
    echo "No replicates detected: Running single-sample workflow for ${SAMPLE_NAME}"
    
    process_replicate "$SAMPLE_NAME"
    
    echo "Skipping Step 4b (Replicate Comparison)."
    echo "Step 5 (single): Generating .shape files from single profile..."
    Rscript "${SHAPE_SCRIPT_DIR}/05_average_reactivities.R" --infiles "${SAMPLE_OUT_DIR}/${SAMPLE_NAME}_profile_annotated.txt" --outdir "$SHAPE_DIR"
fi


# --- Step 6: SHAPE-informed RNAfold ---
echo "Step 6: Running SHAPE-informed RNAfold for ${SAMPLE_NAME}..."
bash "${SHAPE_SCRIPT_DIR}/06_SHAPEinformed_RNAfold.sh" -f "$REF_RNAFOLD_FASTA" -s "$SHAPE_DIR" -o "$RNAFOLD_OUT_DIR"

# --- Step 7: Extract RNAfold Summaries ---
echo "Step 7a: Extracting pairing probabilities..."
Rscript "${RNASTRUCTURE_DIR}/extract_pairing_prob.R" --indir "$RNAFOLD_OUT_DIR" --outfile "${SAMPLE_OUT_DIR}/${SAMPLE_NAME}_pairing_probabilities.csv" --max-len 130 # Adjust if needed

echo "Step 7b: Extracting fold summaries..."
Rscript "${RNASTRUCTURE_DIR}/fold_summary.R" --indir "$RNAFOLD_OUT_DIR" --outfile "${SAMPLE_OUT_DIR}/${SAMPLE_NAME}_rnafold_summary.csv"

echo "================================================================="
echo "PIPELINE COMPLETE for SAMPLE: ${SAMPLE_NAME}"
echo "NOTE: Step 5b (Cross-Sample Comparison) must be run manually after"
echo "all array jobs for the desired samples are complete."
echo "================================================================="