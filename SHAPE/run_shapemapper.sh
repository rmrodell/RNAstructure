#!/bin/bash

# =========================================================================
# STEP 1: ARGUMENT PARSING AND VALIDATION
# This script expects 5 command-line arguments.
# =========================================================================

# Check if the correct number of arguments was provided
if [ "$#" -ne 5 ]; then
    echo "ERROR: Incorrect number of arguments provided."
    echo "Usage: $0 <path_to_modified.fq> <path_to_untreated.fq> <sample_base_name> <base_output_dir> <ref-fasta-dir>"
    exit 1
fi

# Assign command-line arguments to variables for clarity
MODIFIED_FQ="$1"
UNTREATED_FQ="$2"
SAMPLE_BASE_NAME="$3"
BASE_OUTPUT_DIR="$4"
FASTA_DIR="$5"

# =========================================================================
# STEP 2: STATIC CONFIGURATION
# These variables are less likely to change between runs.
# =========================================================================

# --- Path to the SHAPE-Mapper executable ---
SHAPEMAPPER_EXEC="/oak/stanford/groups/nicolemm/shapemapper2-2.3/shapemapper"

# --- Input FASTA Directory (contains all the .fasta chunk files) ---
# FASTA_DIR="/home/groups/nicolemm/rodell/pool1/fasta_chunks"

# --- Resource Allocation (dynamically set by Slurm's --cpus-per-task) ---
N_PROCS=$(nproc)

# =========================================================================
# STEP 3: EXECUTION LOGIC
# No edits should be needed below this line.
# The script will loop through all .fasta files and run SHAPE-Mapper.
# =========================================================================

echo "Starting SHAPE-Mapper batch processing."
echo "Running with the following parameters:"
echo "  - Modified FQ: ${MODIFIED_FQ}"
echo "  - Untreated FQ: ${UNTREATED_FQ}"
echo "  - Sample Base Name: ${SAMPLE_BASE_NAME}"
echo "  - Input FASTA Dir: ${FASTA_DIR}"
echo "  - Base Output Dir: ${BASE_OUTPUT_DIR}"
echo "  - CPUs to use: ${N_PROCS}"
echo "----------------------------------------------------"

# Create the main logs directory once
LOG_DIR="${BASE_OUTPUT_DIR}/logs"
mkdir -p "${LOG_DIR}"

# Find all .fasta files in the FASTA_DIR and loop through them
for REF_FASTA in "${FASTA_DIR}"/*.fasta; do

    # --- Dynamic Name Generation ---
    # Extract the core chunk name (e.g., "pool1_chunk_1") from the full path
    CHUNK_FILENAME=$(basename "${REF_FASTA}" .fasta)
    # Extract just the chunk identifier (e.g., "chunk_1")
    CHUNK_ID=$(echo "${CHUNK_FILENAME}" | awk -F'_' '{print $(NF-1)"_"$(NF)}')
    # Create a unique sample name for this specific chunk
    SAMPLE_NAME="${SAMPLE_BASE_NAME}_${CHUNK_ID}" # e.g., "Rep1_NoPus_chunk_1"
    # Define unique output and log paths for this chunk
    OUTPUT_DIR="${BASE_OUTPUT_DIR}/${SAMPLE_NAME}"
    LOG_FILE="${LOG_DIR}/${SAMPLE_NAME}_shapemapper.log"

    # --- Run Analysis for the Current Chunk ---
    echo "Processing chunk: ${CHUNK_FILENAME}"
    echo " -> Sample Name: ${SAMPLE_NAME}"
    echo " -> Output will be in: ${OUTPUT_DIR}"
    echo " -> Log file will be: ${LOG_FILE}"

    # Ensure the specific output directory for this chunk exists
    mkdir -p "${OUTPUT_DIR}"

    # Execute the SHAPE-Mapper command with dynamically generated names
    ${SHAPEMAPPER_EXEC} \
        --name "${SAMPLE_NAME}" \
        --target "${REF_FASTA}" \
        --out "${OUTPUT_DIR}" \
        --log "${LOG_FILE}" \
        --modified --U "${MODIFIED_FQ}" \
        --untreated --U "${UNTREATED_FQ}" \
        --overwrite \
        --nproc "${N_PROCS}" \
        --min-depth 1000 \
        --window-to-trim 10 \
        --per-read-histograms \
        --indiv-norm \
        --output-counted-mutations \
        --output-aligned-reads

    echo "Finished processing chunk: ${CHUNK_FILENAME}"
    echo "----------------------------------------------------"

done

echo "All chunks processed. Batch analysis complete."