#!/bin/bash
#
# SLURM job array submitter for the SHAPE-MaP pipeline.
# This script submits a single array job, where each task processes one sample.
#

# Stop the script if any command fails
set -e
set -o pipefail

################################################################################
### USER-CONFIGURABLE VARIABLES                                            ###
################################################################################

# --- SLURM JOB RESOURCE CONFIGURATION (for EACH task in the array) ---
CPU_PER_JOB=16
MEM_PER_JOB="32G"
TIME_LIMIT="24:00:00"
PARTITION="normal"
MAIL_USER="rodell@stanford.edu" # <-- IMPORTANT: SET YOUR EMAIL

# --- PATH CONFIGURATION ---
BAM_SOURCE_DIR="/scratch/groups/nicolemm/rodell/basecalling/InVitro_SHAPE_Rep2_20250731/final_results/"
SAMPLE_MAP_FILE="/scratch/groups/nicolemm/rodell/SHAPE/InVitro_Rep2_20250731/samples.txt"
TOP_LEVEL_OUTPUT_DIR="/scratch/groups/nicolemm/rodell/SHAPE/InVitro_Rep2_20250731"
PIPELINE_SCRIPT_PATH="/scratch/groups/nicolemm/rodell/SHAPE/process_single_file.sh"

################################################################################
### BATCH SUBMISSION LOGIC                                                 ###
################################################################################

echo "--- Preparing SLURM Job Array Submission ---"

# Define and create directories
FINAL_FASTQ_DIR="${TOP_LEVEL_OUTPUT_DIR}/final_deduplicated_fastq"
SLURM_LOGS_DIR="${TOP_LEVEL_OUTPUT_DIR}/slurm_logs"
mkdir -p "$FINAL_FASTQ_DIR"
mkdir -p "$SLURM_LOGS_DIR"

echo "All outputs will be placed within: $TOP_LEVEL_OUTPUT_DIR"
echo "SLURM output/error logs will be saved in: $SLURM_LOGS_DIR"

# Check if the pipeline script is executable
if [ ! -x "$PIPELINE_SCRIPT_PATH" ]; then
    echo "Error: Pipeline script '$PIPELINE_SCRIPT_PATH' not found or not executable."
    exit 1
fi

# Determine the number of jobs for the array
N_JOBS=$(wc -l < "$SAMPLE_MAP_FILE")
if [ "$N_JOBS" -eq 0 ]; then
    echo "Error: Sample map file '$SAMPLE_MAP_FILE' is empty."
    exit 1
fi
echo "Found $N_JOBS samples to process."

# Submit the array job
echo "Submitting array job..."
sbatch \
    --job-name="SHAPE_prep" \
    --array=1-"$N_JOBS" \
    --cpus-per-task="$CPU_PER_JOB" \
    --mem="$MEM_PER_JOB" \
    --time="$TIME_LIMIT" \
    --partition="$PARTITION" \
    --mail-type=BEGIN,END,FAIL \
    --mail-user="$MAIL_USER" \
    --output="${SLURM_LOGS_DIR}/job_%A_task_%a.out" \
    --error="${SLURM_LOGS_DIR}/job_%A_task_%a.err" \
    "$PIPELINE_SCRIPT_PATH" "$SAMPLE_MAP_FILE" "$BAM_SOURCE_DIR" "$TOP_LEVEL_OUTPUT_DIR"

echo "--- Job array submitted successfully! ---"