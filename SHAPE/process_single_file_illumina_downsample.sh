#!/bin/bash
#
# A generalized pipeline for processing a paired-end Illumina SHAPE-MaP sample.
# This script is designed to be called as a SLURM array task but can also be run
# in a standalone mode for testing a single sample.
#
# It performs UMI extraction, adapter trimming, mapping, and deduplication,
# producing a final FASTQ file suitable for input into shape-mapper2.
#

# Stop the script if any command fails
set -e
set -o pipefail

# --- 0. Default Variables ---
# This variable is intentionally empty. It will be populated by the -t flag.
# If it remains empty, the script defaults to TASK_ID=1 in standalone mode.
MANUAL_TASK_ID=""

# --- Usage and Argument Parsing ---
usage() {
    echo "Usage: $0 -m <map.tsv> -o <out_dir> -r <ref_index> [-t <task_id>]"
    echo "       $0 --map_file <map.tsv> --output_dir <out_dir> --ref_index <ref_index> [--task_id <task_id>]"
    echo ""
    echo "  -m, --map_file      <map.tsv>      : Path to the tab-separated file mapping sample IDs to R1/R2 files."
    echo "  -o, --output_dir    <out_dir>      : The root directory where all output will be stored."
    echo "  -r, --ref_index     <ref_index>    : Full path prefix for the Bowtie2 index."
    echo "  -t, --task_id       <task_id>      : (Optional) Manual task ID to run a single sample outside of a SLURM array."
    echo "  -h, --help                         : Display this help message."
    exit 1
}

# Use getopt to parse long and short options
OPTS=$(getopt -o hm:o:r:t: -l help,map_file:,output_dir:,ref_index:,task_id: --name "$0" -- "$@")
if [ $? != 0 ]; then echo "Failed to parse options." >&2; usage; fi
eval set -- "$OPTS"

while true; do
    case "$1" in
        -h | --help ) usage ;;
        -m | --map_file ) SAMPLE_MAP_FILE="$2"; shift 2 ;;
        -o | --output_dir ) TOP_LEVEL_OUTPUT_DIR="$2"; shift 2 ;;
        -r | --ref_index ) BOWTIE2_INDEX="$2"; shift 2 ;;
        -t | --task_id ) MANUAL_TASK_ID="$2"; shift 2 ;;
        -- ) shift; break ;;
        * ) break ;;
    esac
done

# Check for mandatory arguments
if [ -z "$SAMPLE_MAP_FILE" ] || [ -z "$TOP_LEVEL_OUTPUT_DIR" ] || [ -z "$BOWTIE2_INDEX" ]; then
    echo "Error: Missing one or more required arguments."
    usage
fi

# --- 1. SLURM/Standalone Task Setup ---
# Prioritize SLURM array task ID, fall back to manual ID, then default to 1.
# The ${VAR:-default} syntax provides the default value if VAR is empty or unset.
TASK_ID=${SLURM_ARRAY_TASK_ID:-${MANUAL_TASK_ID:-1}}

# Echo the assigned arguments to the log for traceability
echo "--- Script Arguments Received ---"
echo "Sample Map File:          ${SAMPLE_MAP_FILE}"
echo "Top-Level Output Dir:     ${TOP_LEVEL_OUTPUT_DIR}"
echo "Bowtie2 Index:            ${BOWTIE2_INDEX}"
echo "Effective Task ID:        ${TASK_ID}"
echo "---------------------------------"
echo ""

# Read the line from the sample map corresponding to the task ID
TASK_INFO=$(sed -n "${TASK_ID}p" "$SAMPLE_MAP_FILE")
if [ -z "$TASK_INFO" ]; then
    echo "Error: No data found in '$SAMPLE_MAP_FILE' for task ID $TASK_ID."
    exit 1
fi

# Parse the sample ID and file paths (Tab Separated)
IFS=$'\t' read -r SAMPLE_ID RAW_R1 RAW_R2 <<< "$TASK_INFO"

# Validate that paths were found and files exist
if [ -z "$SAMPLE_ID" ] || [ -z "$RAW_R1" ] || [ -z "$RAW_R2" ]; then
    echo "Error: Malformed line in sample map for task ID $TASK_ID. Expected format: SampleID<tab>PathToR1<tab>PathToR2"
    exit 1
fi
if [ ! -f "$RAW_R1" ]; then echo "Error: R1 file not found: $RAW_R1"; exit 1; fi
if [ ! -f "$RAW_R2" ]; then echo "Error: R2 file not found: $RAW_R2"; exit 1; fi


# --- 2. Load Modules & Define Variables ---
ml biology py-cutadapt/1.18_py36 samtools bowtie2/2.3.4.1
ml python/3.6.1
ml devel java/11

POOL_SENSE_ADAPTER_3PRIME="CACTCGGGCACCAAGGAC"
UMI_PATTERN="NNNNNNNNNN"

# SLURM-aware thread count
if [ -n "$SLURM_CPUS_PER_TASK" ]; then THREADS="$SLURM_CPUS_PER_TASK"; else THREADS=4; fi

# --- 3. Directory Setup & Helper Functions ---
PROJECT_DIR="${TOP_LEVEL_OUTPUT_DIR}/${SAMPLE_ID}"
FINAL_COMMON_DIR="${TOP_LEVEL_OUTPUT_DIR}/final_deduplicated_fastq"
TMP_DIR="${PROJECT_DIR}/tmp"
FINAL_DIR="${PROJECT_DIR}/final"
LOGS_DIR="${PROJECT_DIR}/logs"

mkdir -p "$TMP_DIR" "$FINAL_DIR" "$LOGS_DIR" "$FINAL_COMMON_DIR"

# Define log files and redirect all output
PIPELINE_LOG="${LOGS_DIR}/${SAMPLE_ID}.pipeline_run.log"
exec > >(tee -a "${PIPELINE_LOG}") 2>&1

METRICS_FILE="${LOGS_DIR}/${SAMPLE_ID}.run_metrics.tsv"
echo -e "Step\tSample\tFile\tRecord_Count\tFile_Size" > "$METRICS_FILE"

log_message() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [Task:${TASK_ID}] [$SAMPLE_ID] $1"; }
format_duration() { local s=$1; echo "$((s/60))m $((s%60))s"; }
track_metrics() {
    local step_name=$1; local file_path=$2; local count; local size;
    if [ ! -s "$file_path" ]; then
        echo -e "${step_name}\t${SAMPLE_ID}\t$(basename "$file_path")\t0\t0B" >> "$METRICS_FILE"; return
    fi
    size=$(du -h "$file_path" | cut -f1)
    if [[ "$file_path" == *.fastq.gz || "$file_path" == *.fq.gz ]]; then
        count=$(echo "$(zcat "$file_path" | wc -l) / 4" | bc)
    elif [[ "$file_path" == *.bam || "$file_path" == *.sam ]]; then
        count=$(samtools view -@ "$THREADS" -c "$file_path")
    else
        count="N/A"
    fi
    echo -e "${step_name}\t${SAMPLE_ID}\t$(basename "$file_path")\t${count}\t${size}" >> "$METRICS_FILE"
}

# --- 4. Define Intermediate File Paths ---
R1_UMI="${TMP_DIR}/${SAMPLE_ID}_R1_umi.fastq.gz"
R2_UMI="${TMP_DIR}/${SAMPLE_ID}_R2_umi.fastq.gz"
TRIMMED_R1="${TMP_DIR}/${SAMPLE_ID}_R1_trimmed.fastq.gz"
MAPPED_SAM="${TMP_DIR}/${SAMPLE_ID}_mapped.sam"
MAPPED_BAM="${TMP_DIR}/${SAMPLE_ID}_mapped_sorted.bam"
DEDUP_BAM="${FINAL_DIR}/${SAMPLE_ID}_deduplicated.bam"
DEDUP_FASTQ="${FINAL_DIR}/${SAMPLE_ID}_deduplicated.fastq"

################################################################################
### PIPELINE STEPS                                                           ###
################################################################################

log_message "--- Starting Illumina SHAPE-MaP pipeline for ${SAMPLE_ID} ---"
log_message "--- Using ${THREADS} threads for parallel tasks ---"
log_message "--- Input R1: ${RAW_R1} ---"
log_message "--- Input R2: ${RAW_R2} ---"

# --- Step 1: Extract UMIs from R2 and add to R1 headers ---
log_message "Step 1: Extracting UMIs from R2 to R1..."
start_time=$(date +%s)
umi_tools extract \
    -I "$RAW_R2" \
    --read2-in="$RAW_R1" \
    --bc-pattern="$UMI_PATTERN" \
    --extract-method=string \
    --read2-out="$R1_UMI" \
    -L "${LOGS_DIR}/${SAMPLE_ID}_umi_extract.log" \
    | gzip > "$R2_UMI"
track_metrics "1_R1_with_UMI" "$R1_UMI"
end_time=$(date +%s); log_message "Step 1 finished. Duration: $(format_duration $((end_time - start_time)))"

# --- Step 2: Trim 3' adapter from R1 reads ---
log_message "Step 2: Trimming 3' adapter from R1 reads..."
start_time=$(date +%s)
cutadapt \
    -m 125 \
    -O 12 \
    --cores="$THREADS" \
    -a "$POOL_SENSE_ADAPTER_3PRIME" \
    -o "$TRIMMED_R1" \
    "$R1_UMI" \
    > "${LOGS_DIR}/${SAMPLE_ID}_cutadapt_trim.log"
track_metrics "2_Trimmed_R1" "$TRIMMED_R1"
end_time=$(date +%s); log_message "Step 2 finished. Duration: $(format_duration $((end_time - start_time)))"

# --- Step 3: Map reads with Bowtie2 ---
log_message "Step 3: Mapping reads with Bowtie2..."
start_time=$(date +%s)
bowtie2 \
    -x "$BOWTIE2_INDEX" \
    -U "$TRIMMED_R1" \
    -S "$MAPPED_SAM" \
    -p "$THREADS" \
    --local \
    --sensitive-local \
    --ignore-quals \
    --no-unal \
    --mp 3,1 \
    --rdg 5,1 \
    --rfg 5,1 \
    --dpad 30 \
    2> "${LOGS_DIR}/${SAMPLE_ID}_bowtie2_mapping.log"
track_metrics "3_Mapped_SAM" "$MAPPED_SAM"
end_time=$(date +%s); log_message "Step 3 finished. Duration: $(format_duration $((end_time - start_time)))"

# --- Step 4: Convert SAM to sorted BAM and index ---
log_message "Step 4: Converting SAM to sorted BAM..."
start_time=$(date +%s)
samtools view -@ "$THREADS" -b "$MAPPED_SAM" | samtools sort -@ "$THREADS" -o "$MAPPED_BAM"
samtools index "$MAPPED_BAM"
track_metrics "4_Sorted_BAM" "$MAPPED_BAM"
end_time=$(date +%s); log_message "Step 4 finished. Duration: $(format_duration $((end_time - start_time)))"

# --- Step 5: Downsample BAM to 20M reads ---
log_message "Step 5: Downsampling mapped BAM to a maximum of 20M reads..."
start_time=$(date +%s)
DOWNSAMPLED_BAM="${TMP_DIR}/${SAMPLE_ID}_mapped_sorted_downsampled.bam"
TARGET_READS=20000000
TOTAL_READS=$(samtools view -c "$MAPPED_BAM")

if [ "$TOTAL_READS" -le "$TARGET_READS" ]; then
    log_message "  - Total reads ($TOTAL_READS) is less than or equal to target ($TARGET_READS). Skipping downsampling."
    cp "$MAPPED_BAM" "$DOWNSAMPLED_BAM"
else
    # Calculate fraction needed for samtools view -s INT.FRAC
    # Using bc for floating point arithmetic
    FRACTION=$(echo "scale=8; $TARGET_READS / $TOTAL_READS" | bc)
    log_message "  - Downsampling $TOTAL_READS reads to ~$TARGET_READS reads (fraction: $FRACTION)."
    # The -s seed is set using the task ID for reproducibility. The "##.}" removes the leading "0." from the fraction.
    samtools view -@ "$THREADS" -b -s "${TASK_ID}.${FRACTION##*.}" "$MAPPED_BAM" > "$DOWNSAMPLED_BAM"
    samtools index "$DOWNSAMPLED_BAM"
fi
track_metrics "5_Downsampled_BAM" "$DOWNSAMPLED_BAM"
end_time=$(date +%s); log_message "Step 5 finished. Duration: $(format_duration $((end_time - start_time)))"


# --- Step 6: Deduplicate reads with UMI-tools ---
log_message "Step 6: Deduplicating reads with UMI-tools and indexing..."
start_time=$(date +%s)
$HOME/UMICollapse/umicollapse bam \
    -i "$MAPPED_BAM" \
    -o "$DEDUP_BAM"
# umi_tools dedup \
#     --method directional \
#     -I "$DOWNSAMPLED_BAM" \
#     -S "$DEDUP_BAM" \
#     -L "${LOGS_DIR}/${SAMPLE_ID}_umi_tools_dedup.log"
samtools index "$DEDUP_BAM"
track_metrics "6_Deduplicated_BAM" "$DEDUP_BAM"
end_time=$(date +%s); log_message "Step 6 finished. Duration: $(format_duration $((end_time - start_time)))"

# --- Step 7: Convert deduplicated BAM to final FASTQ for SHAPE-Mapper ---
log_message "Step 7: Converting deduplicated BAM to final FASTQ..."
start_time=$(date +%s)
samtools bam2fq -@ "$THREADS" "$DEDUP_BAM" > "$DEDUP_FASTQ"
track_metrics "7_Final_Deduplicated_FASTQ" "$DEDUP_FASTQ"
end_time=$(date +%s); log_message "Step 7 finished. Duration: $(format_duration $((end_time - start_time)))"

# --- Step 8: Copy final FASTQ to common output directory and verify ---
log_message "Step 8: Copying and verifying final FASTQ..."
start_time=$(date +%s)
DEST_FASTQ_PATH="${FINAL_COMMON_DIR}/${SAMPLE_ID}_deduplicated.fastq"
cp "$DEDUP_FASTQ" "$DEST_FASTQ_PATH"
source_md5=$(md5sum "$DEDUP_FASTQ" | awk '{print $1}')
dest_md5=$(md5sum "$DEST_FASTQ_PATH" | awk '{print $1}')
if [ "$source_md5" == "$dest_md5" ]; then
    log_message "  - SUCCESS: Checksums match. File integrity confirmed."
else
    log_message "  - FATAL ERROR: Checksums DO NOT MATCH!"
    exit 1
fi
end_time=$(date +%s); log_message "Step 8 finished. Duration: $(format_duration $((end_time - start_time)))"

# --- Step 9: Clean up intermediate files ---
log_message "Step 9: Cleaning up temporary directory..."
rm -rf "$TMP_DIR"

# --- Step 10: Display Final Metrics Report ---
log_message "Step 10: Displaying run metrics summary for ${SAMPLE_ID}:"
column -t -s $'\t' "$METRICS_FILE"

log_message "--- Pipeline for ${SAMPLE_ID} finished successfully! ---"