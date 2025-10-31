#!/bin/bash

# =========================================================================
# Description:
#   This script serves as a wrapper to run the SHAPE-Mapper 2 executable in a batch
#   mode. It iterates through all reference FASTA files (.fa or .fasta) in a
#   specified directory, running a separate analysis for each one.
#
# shapemapper2-2.3 must be installed into an executeable path; update below if needed
#
# Arguments:
#   -m, --modified   <path>   Path to the modified sample FASTQ file. (Required)
#   -u, --untreated  <path>   Path to the untreated sample FASTQ file. (Required)
#   -r, --ref-fasta  <path>   Path to the directory containing reference .fa/.fasta files. (Required)
#   -n, --name       <string> Base name for the sample series, used for output files. (Required)
#   -o, --output-dir <path>   Path to the base directory for all outputs. (Required)
#   -h, --help                Display this help message and exit.
#
# FASTA File Naming Convention:
#   The script assumes reference files are named like 'pool1_chunk_1.fasta',
#   'pool1_chunk_2.fasta', etc. It extracts the 'chunk_#' part to create
#   unique output names.
#
# Example Usage:
#   ./run_shapemapper_batch.sh \
#       -m /path/to/my/modified.fastq.gz \
#       -u /path/to/my/untreated.fastq.gz \
#       -r /path/to/my/fasta_chunks_dir \
#       -n MySample_Rep1 \
#       -o /path/to/output_directory
#
# ===================================================================================

# --- Define a function for displaying usage information ---
usage() {
    echo "Usage: $0 -m <modified.fq> -u <untreated.fq> -r <ref_fasta_dir> -n <sample_name> -o <output_dir>"
    echo ""
    echo "Required Arguments:"
    echo "  -m, --modified   <path>   Path to the modified FASTQ file."
    echo "  -u, --untreated  <path>   Path to the untreated FASTQ file."
    echo "  -r, --ref-fasta  <path>   Directory containing reference FASTA files (.fa, .fasta)."
    echo "  -n, --name       <string> Base name for the sample series."
    echo "  -o, --output-dir <path>   Base directory for all outputs and the main log file."
    echo ""
    echo "Optional Arguments:"
    echo "  -h, --help                Display this help message and exit."
    exit 1
}

# --- Initialize variables ---
MODIFIED_FQ=""
UNTREATED_FQ=""
FASTA_DIR=""
SAMPLE_BASE_NAME=""
BASE_OUTPUT_DIR=""

# =========================================================================
# STEP 1: ARGUMENT PARSING AND VALIDATION
# =========================================================================

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -m|--modified) MODIFIED_FQ="$2"; shift ;;
        -u|--untreated) UNTREATED_FQ="$2"; shift ;;
        -r|--ref-fasta) FASTA_DIR="$2"; shift ;;
        -n|--name) SAMPLE_BASE_NAME="$2"; shift ;;
        -o|--output-dir) BASE_OUTPUT_DIR="$2"; shift ;;
        -h|--help) usage ;;
        *) echo "ERROR: Unknown parameter passed: $1"; usage ;;
    esac
    shift
done

# Check for missing required arguments
if [[ -z "$MODIFIED_FQ" || -z "$UNTREATED_FQ" || -z "$FASTA_DIR" || -z "$SAMPLE_BASE_NAME" || -z "$BASE_OUTPUT_DIR" ]]; then
    echo "ERROR: One or more required arguments are missing."
    usage
fi

# Check if the reference FASTA directory exists
if [[ ! -d "$FASTA_DIR" ]]; then
    echo "ERROR: Reference FASTA directory not found at: ${FASTA_DIR}"
    exit 1
fi

# Create the main output and logs directories
LOG_DIR="${BASE_OUTPUT_DIR}/logs"
mkdir -p "${BASE_OUTPUT_DIR}"
mkdir -p "${LOG_DIR}"

LOG_FILE="${LOG_DIR}/${SAMPLE_BASE_NAME}_batch_run.log"

# Redirect all stdout and stderr from this point forward ONLY to the log file
exec &> "${LOG_FILE}"

# =========================================================================
# STEP 2: STATIC CONFIGURATION
# =========================================================================

# --- Path to the SHAPE-Mapper executable ---
SHAPEMAPPER_EXEC="/scratch/groups/nicolemm/rodell/shapemapper2-2.3/shapemapper"

# --- Resource Allocation (dynamically set by Slurm's --cpus-per-task) ---
N_PROCS=$(nproc)

echo "===================================================="
echo "Starting SHAPE-Mapper BATCH analysis on $(date)"
echo "----------------------------------------------------"
echo "Parameters provided:"
echo "  - Modified FASTQ:     ${MODIFIED_FQ}"
echo "  - Untreated FASTQ:    ${UNTREATED_FQ}"
echo "  - Reference FASTA Dir:${FASTA_DIR}"
echo "  - Sample Base Name:   ${SAMPLE_BASE_NAME}"
echo "  - Base Output Dir:    ${BASE_OUTPUT_DIR}"
echo ""
echo "Configuration:"
echo "  - Log File:           ${LOG_FILE}"
echo "  - SHAPE-Mapper Exec:  ${SHAPEMAPPER_EXEC}"
echo "  - CPUs to use:        ${N_PROCS}"
echo "===================================================="

# =========================================================================
# STEP 3: EXECUTION LOGIC
# No edits should be needed below this line.
# The script will loop through all .fasta files and run SHAPE-Mapper.
# =========================================================================

# # Find all .fasta files in the FASTA_DIR and loop through them
# for REF_FASTA in "${FASTA_DIR}"/*.fasta; do

#     # --- Dynamic Name Generation ---
#     # Extract the core chunk name (e.g., "pool1_chunk_1") from the full path
#     CHUNK_FILENAME=$(basename "${REF_FASTA}" .fasta)
#     # Extract just the chunk identifier (e.g., "chunk_1")
#     CHUNK_ID=$(echo "${CHUNK_FILENAME}" | cut -d'_' -f2,3)
#     # Create a unique sample name for this specific chunk
#     SAMPLE_NAME="${SAMPLE_BASE_NAME}_${CHUNK_ID}" # e.g., "Rep1_NoPus_chunk_1"
#     # Define unique output and log paths for this chunk
#     OUTPUT_DIR="${BASE_OUTPUT_DIR}/${SAMPLE_NAME}"
#     LOG_FILE="${LOG_DIR}/${SAMPLE_NAME}_shapemapper.log"

#     # --- Run Analysis for the Current Chunk ---
#     echo "Processing chunk: ${CHUNK_FILENAME}"
#     echo " -> Sample Name: ${SAMPLE_NAME}"
#     echo " -> Output will be in: ${OUTPUT_DIR}"
#     echo " -> Log file will be: ${LOG_FILE}"

#     # Ensure the specific output directory for this chunk exists
#     mkdir -p "${OUTPUT_DIR}"

#     # Execute the SHAPE-Mapper command with dynamically generated names
#     ${SHAPEMAPPER_EXEC} \
#         --name "${SAMPLE_NAME}" \
#         --target "${REF_FASTA}" \
#         --out "${OUTPUT_DIR}" \
#         --log "${LOG_FILE}" \
#         --modified --U "${MODIFIED_FQ}" \
#         --untreated --U "${UNTREATED_FQ}" \
#         --overwrite \
#         --nproc "${N_PROCS}" \
#         --min-depth 1000 \
#         --window-to-trim 10 \
#         --per-read-histograms \
#         --indiv-norm \
#         --output-counted-mutations \
#         --output-aligned-reads

#     echo "Finished processing chunk: ${CHUNK_FILENAME}"
#     echo "----------------------------------------------------"

# done


# Find all .fa and .fasta files and loop through them
shopt -s nullglob
FASTA_FILES=("${FASTA_DIR}"/*.{fa,fasta})
shopt -u nullglob

if [ ${#FASTA_FILES[@]} -eq 0 ]; then
    echo "ERROR: No .fa or .fasta files found in ${FASTA_DIR}. Exiting."
    exit 1
fi

echo "Found ${#FASTA_FILES[@]} reference files to process."

# Loop through all .fa / .fasta files
for REF_FASTA in "${FASTA_FILES[@]}"; do

    # --- Dynamic Name Generation ---
    FILENAME=$(basename "${REF_FASTA}")
    CHUNK_FILENAME="${FILENAME%.*}" # Removes extension (.fa or .fasta)
    CHUNK_ID=$(echo "${CHUNK_FILENAME}" | cut -d'_' -f2,3) # Assumes 'pool_chunk_id' format
    SAMPLE_NAME="${SAMPLE_BASE_NAME}_${CHUNK_ID}"
    OUTPUT_SUBDIR="${BASE_OUTPUT_DIR}/${SAMPLE_NAME}"
    CHUNK_LOG_FILE="${OUTPUT_SUBDIR}/${SAMPLE_NAME}_shapemapper.log"

    # --- Run Analysis for the Current Chunk ---
    echo "----------------------------------------------------"
    echo "Processing chunk: ${FILENAME}"
    echo " -> Sample Name: ${SAMPLE_NAME}"
    echo " -> Output will be in: ${OUTPUT_SUBDIR}"

    mkdir -p "${OUTPUT_SUBDIR}"

    # Execute the SHAPE-Mapper command with dynamically generated names
    ${SHAPEMAPPER_EXEC} \
        --name "${SAMPLE_NAME}" \
        --target "${REF_FASTA}" \
        --out "${OUTPUT_SUBDIR}" \
        --modified --U "${MODIFIED_FQ}" \
        --untreated --U "${UNTREATED_FQ}" \
        --log "${CHUNK_LOG_FILE}" \
        --overwrite \
        --nproc "${N_PROCS}" \
        --min-depth 1000 \
        --window-to-trim 10 \
        --per-read-histograms \
        --indiv-norm \
        --output-counted-mutations \
        --output-aligned-reads

    echo "Finished processing chunk: ${FILENAME}"
    echo "----------------------------------------------------"

done

echo "All chunks processed. Batch analysis complete."
echo "Full log is available at: ${LOG_FILE}"