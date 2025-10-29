#!/bin/bash
#
# A resume script to complete the pipeline for a SINGLE SHAPE-MaP sample.
# This script starts from the deduplication step (Step 8).
#

# --- 0. Validate Input Arguments ---
if [ "$#" -ne 3 ]; then
    echo "Error: Incorrect number of arguments provided."
    echo ""
    echo "Usage: $0 <sample_id> <project_dir> <final_common_dir>"
    echo ""
    echo "  <sample_id>        : The unique identifier for the sample (e.g., HEK293T_DMSO_Rep1)."
    echo "  <project_dir>      : The sample-specific output directory containing the 'tmp' folder."
    echo "  <final_common_dir> : The root directory where the final FASTQ will be placed."
    exit 1
fi

# Assign arguments to variables
SAMPLE_ID="$1"
PROJECT_DIR="$2"
FINAL_COMMON_DIR="$3"


echo "--- Resume Script Initialized ---"
echo "Sample ID:                ${SAMPLE_ID}"
echo "Project Directory:        ${PROJECT_DIR}"
echo "Final FASTQ Directory:    ${FINAL_COMMON_DIR}"
echo "---------------------------------"
echo ""

# Load modules
ml biology py-cutadapt/1.18_py36 samtools bowtie2/2.3.4.1
ml python/3.6.1

# Stop the script if any command fails
set -e
set -o pipefail

################################################################################
###             HELPER FUNCTIONS AND DIRECTORY SETUP                         ###
################################################################################

# SLURM-aware thread count
if [ -n "$SLURM_CPUS_PER_TASK" ]; then THREADS="$SLURM_CPUS_PER_TASK"; else THREADS=1; fi

# Use SLURM_JOB_ID for logging consistency
log_message() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [Job:${SLURM_JOB_ID}] [${SAMPLE_ID}] $1"; }
format_duration() { local s=$1; echo "$((s/60))m $((s%60))s"; }

# --- Directory and Log File Definitions (must match original) ---
TMP_DIR="${PROJECT_DIR}/tmp"
FINAL_DIR="${PROJECT_DIR}/final"
LOGS_DIR="${PROJECT_DIR}/logs"
REPORTS_DIR="${PROJECT_DIR}/reports"

# Ensure directories exist
mkdir -p "$TMP_DIR" "$FINAL_DIR" "$LOGS_DIR" "$REPORTS_DIR" "$FINAL_COMMON_DIR"

METRICS_FILE="${REPORTS_DIR}/${SAMPLE_ID}.run_metrics.tsv"
PIPELINE_LOG="${LOGS_DIR}/${SAMPLE_ID}.pipeline_run.log"

# Append to existing logs instead of overwriting
exec >> >(tee -a "${PIPELINE_LOG}") 2>&1

# Define intermediate and final file paths (must match original)
MAPPED_BAM="${TMP_DIR}/${SAMPLE_ID}_mapped_sorted.bam"
DEDUP_BAM="${FINAL_DIR}/${SAMPLE_ID}_deduplicated.bam"
DEDUP_FASTQ="${FINAL_DIR}/${SAMPLE_ID}_deduplicated_for_shapemapper.fastq"
DEDUP_LOG="${LOGS_DIR}/${SAMPLE_ID}_umi_tools_dedup.log"


################################################################################
### PIPELINE RESUMPTION                                                      ###
################################################################################

log_message "--- Resuming pipeline from Step 8 for ${SAMPLE_ID} ---"
log_message "--- Using ${THREADS} threads for parallel tasks ---"

# --- VALIDATE PREREQUISITE FILE ---
if [ ! -f "$MAPPED_BAM" ]; then
    log_message "FATAL ERROR: Prerequisite file for Step 8 not found!"
    log_message "Expected to find: ${MAPPED_BAM}"
    log_message "Cannot resume pipeline. Please check the path and the output of the previous run."
    exit 1
fi
log_message "Prerequisite file found: ${MAPPED_BAM}"
log_message "Skipping Steps 1-7."

# --- 8. Deduplicate Reads with UMI-tools ---
log_message "Step 8: Deduplicating reads with UMI-tools and indexing..."
start_time=$(date +%s)
umi_tools dedup \
    --method directional \
    -I "$MAPPED_BAM" \
    -S "$DEDUP_BAM" \
    -L "$DEDUP_LOG"
samtools index "$DEDUP_BAM"
# The original script's track_metrics is omitted here, as it's not essential for resumption.
# You could add it back if you want to regenerate the metrics for this step.
end_time=$(date +%s); log_message "Step 8 finished. Duration: $(format_duration $((end_time - start_time)))"

# --- 9. Convert Deduplicated BAM to FASTQ for SHAPE-Mapper ---
log_message "Step 9: Converting deduplicated BAM to final FASTQ..."
start_time=$(date +%s)
samtools bam2fq -@ "$THREADS" "$DEDUP_BAM" > "$DEDUP_FASTQ"
end_time=$(date +%s); log_message "Step 9 finished. Duration: $(format_duration $((end_time - start_time)))"

# --- Step 10: Copy Final FASTQ to Common Output Directory ---
log_message "Step 10: Copying and verifying final FASTQ..."
start_time=$(date +%s)

DEST_FASTQ_PATH="${FINAL_COMMON_DIR}/${SAMPLE_ID}.final.fastq"

cp "$DEDUP_FASTQ" "$DEST_FASTQ_PATH"
log_message "  - Copied to: ${DEST_FASTQ_PATH}"

log_message "  - Verifying file integrity with md5sum..."
source_md5=$(md5sum "$DEDUP_FASTQ" | awk '{print $1}')
dest_md5=$(md5sum "$DEST_FASTQ_PATH" | awk '{print $1}')

if [ "$source_md5" == "$dest_md5" ]; then
    log_message "  - SUCCESS: Checksums match. File integrity confirmed."
else
    log_message "  - FATAL ERROR: Checksums DO NOT MATCH!"
    log_message "    - Source:      ${source_md5}"
    log_message "    - Destination: ${dest_md5}"
    log_message "  - Halting pipeline to prevent use of corrupted data."
    exit 1
fi
end_time=$(date +%s); log_message "Step 10 finished. Duration: $(format_duration $((end_time - start_time)))"

# --- Step 11: Clean up intermediate files ---
log_message "Step 11: Cleaning up temporary directory..."
rm -rf "$TMP_DIR"
log_message "  - Removed: ${TMP_DIR}"

log_message "--- Pipeline for ${SAMPLE_ID} finished successfully! ---"

# --- Step 12: (Optional) Display Final Metrics Report ---
# Note: This will show an incomplete report as we didn't re-run track_metrics.
# It's left here for consistency but you can ignore its output.
if [ -f "$METRICS_FILE" ]; then
    log_message "Step 12: Displaying run metrics summary for ${SAMPLE_ID}:"
    column -t -s $'\t' "$METRICS_FILE"
fi