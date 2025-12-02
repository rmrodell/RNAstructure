#!/bin/bash
#
# shapemapper_single_target.sh
#
# Formalized script to run a single target sequence through the modular ShapeMapper2 pipeline.
#
# This script automates the following steps:
# 1. Filters SAM files to retain only alignments to a specific target sequence.
# 2. Parses mutations from the filtered SAM files for modified and untreated samples.
# 3. Counts mutations for each sample (now correctly sized for the target).
# 4. Calculates reactivity profiles for the specific target sequence.
# 5. Normalizes the final reactivity profile.
#
# It creates all necessary output directories and captures all command outputs in a log file.

set -e
set -o pipefail

# --- Usage function ---
usage() {
    echo "Usage: $0 -d <shapemapper_dir> -f <fasta_file> -r <target_rna_name> -m <MODIFIED_BAM> -u <UNTREATED_BAM> -o <output_dir> -p <prefix>"
    echo "  -d  Path to the main shapemapper2-2.3 directory."
    echo "  -f  Path to the multi-sequence FASTA reference file."
    echo "  -r  Name of the target RNA sequence (must match FASTA header)."
    echo "  -m  Path to the 'modified' indexed sample BAM file."
    echo "  -u  Path to the 'untreated' indexed (e.g., DMSO) sample BAM file."
    echo "  -o  Path to the main output directory for all results."
    echo "  -p  A prefix for output file names (e.g., 'Pool2_Rep1')."
    exit 1
}

# --- Parse Command-Line Arguments ---
while getopts ":d:f:r:m:u:o:p:" opt; do
    case ${opt} in
        d) SHAPEMAPPER_DIR=$OPTARG ;;
        f) FASTA_FILE=$OPTARG ;;
        r) TARGET_RNA=$OPTARG ;;
        m) MODIFIED_BAM=$OPTARG ;;
        u) UNTREATED_BAM=$OPTARG ;;
        o) OUTPUT_DIR=$OPTARG ;;
        p) PREFIX=$OPTARG ;;
        \?) echo "Invalid option: -$OPTARG" >&2; usage ;;
        :) echo "Option -$OPTARG requires an argument." >&2; usage ;;
    esac
done

# --- Check for completeness of arguments ---
if [ -z "$SHAPEMAPPER_DIR" ] || [ -z "$FASTA_FILE" ] || [ -z "$TARGET_RNA" ] || [ -z "$MODIFIED_BAM" ] || [ -z "$UNTREATED_BAM" ] || [ -z "$OUTPUT_DIR" ] || [ -z "$PREFIX" ]; then
    echo "Error: Missing one or more required arguments."
    usage
fi

# --- Setup Directories and Logging ---
mkdir -p "$OUTPUT_DIR"
LOG_FILE="${OUTPUT_DIR}/${PREFIX}_${TARGET_RNA}_pipeline.log"
exec > >(tee -a "$LOG_FILE") 2>&1

# --- Print Variables for Inspection ---
echo "================================================="
echo "Starting ShapeMapper Pipeline for a Single Target"
echo "================================================="
echo "Date: $(date)"
echo "ShapeMapper2 Directory: ${SHAPEMAPPER_DIR}"
echo "FASTA File: ${FASTA_FILE}"
echo "Target RNA: ${TARGET_RNA}"
echo "Modified BAM: ${MODIFIED_BAM}"
echo "Untreated BAM: ${UNTREATED_BAM}"
echo "Output Directory: ${OUTPUT_DIR}"
echo "File Prefix: ${PREFIX}"
echo "Log File: ${LOG_FILE}"
echo "-------------------------------------------------"


# --- Define Variables for tools and intermediate files ---
SHAPEMAPPER_BIN="${SHAPEMAPPER_DIR}/internals/bin"

# Dynamic file names
MODIFIED_SAM_FILTERED="${OUTPUT_DIR}/${PREFIX}_modified_${TARGET_RNA}.sam"
UNTREATED_SAM_FILTERED="${OUTPUT_DIR}/${PREFIX}_untreated_${TARGET_RNA}.sam"
MODIFIED_MUT="${OUTPUT_DIR}/${PREFIX}_modified_${TARGET_RNA}.mut"
UNTREATED_MUT="${OUTPUT_DIR}/${PREFIX}_untreated_${TARGET_RNA}.mut"
MODIFIED_COUNTS="${OUTPUT_DIR}/${PREFIX}_modified_${TARGET_RNA}_counts.txt"
UNTREATED_COUNTS="${OUTPUT_DIR}/${PREFIX}_untreated_${TARGET_RNA}_counts.txt"
PROFILE="${OUTPUT_DIR}/${PREFIX}_${TARGET_RNA}_profile.txt"

# --- Pipeline Steps ---
ml biology samtools

# 1. Filter BAM files for the target RNA
# This is the CRITICAL step to fix the length mismatch error.
echo "[Step 1] Filtering BAM files for target: ${TARGET_RNA}"
samtools view -h -@ 4 "$MODIFIED_BAM" "$TARGET_RNA" > "$MODIFIED_SAM_FILTERED"
samtools view -h -@ 4 "$UNTREATED_BAM" "$TARGET_RNA" > "$UNTREATED_SAM_FILTERED"
echo "Filtered BAM files created."

# 2. Parse mutations from filtered SAM files
echo "[Step 2] Parsing mutations from filtered SAM files..."
TEMP_PARSER_LOG_MOD="${OUTPUT_DIR}/temp_parser_mod.log"
TEMP_PARSER_LOG_UN="${OUTPUT_DIR}/temp_parser_un.log"

source $GROUP_SCRATCH/rodell/shapemapper2-2.3/internals/paths/bin_paths.sh

echo "  - Parsing modified sample..."
{
    "${SHAPEMAPPER_BIN}/shapemapper_mutation_parser" -i "$MODIFIED_SAM_FILTERED" -o "$MODIFIED_MUT" --input_is_unpaired
} > "$TEMP_PARSER_LOG_MOD" 2>&1
rm "$TEMP_PARSER_LOG_MOD" # Delete the temporary verbose log

echo "  - Parsing untreated sample..."
{
    "${SHAPEMAPPER_BIN}/shapemapper_mutation_parser" -i "$UNTREATED_SAM_FILTERED" -o "$UNTREATED_MUT" --input_is_unpaired
} > "$TEMP_PARSER_LOG_UN" 2>&1
rm "$TEMP_PARSER_LOG_UN" # Delete the temporary verbose log

echo "Mutation parsing complete. Temporary logs deleted."

# 3. Count mutations
echo "[Step 3] Counting mutations..."
"${SHAPEMAPPER_BIN}/shapemapper_mutation_counter" -i "$MODIFIED_MUT" -c "$MODIFIED_COUNTS"
"${SHAPEMAPPER_BIN}/shapemapper_mutation_counter" -i "$UNTREATED_MUT" -c "$UNTREATED_COUNTS"
echo "Mutation counting complete."

# 4. Calculate reactivity profile
echo "[Step 4] Calculating reactivity profile..."
"${SHAPEMAPPER_BIN}/make_reactivity_profiles.py" \
    --fa "$FASTA_FILE" \
    --rna "$TARGET_RNA" \
    --counts "$MODIFIED_COUNTS" "$UNTREATED_COUNTS" \
    --out "$PROFILE"
echo "Raw profile created at: ${PROFILE}"

# 5. Normalize profile
echo "[Step 5] Normalizing reactivity profile..."
"${SHAPEMAPPER_BIN}/normalize_profiles.py" --tonorm "$PROFILE"
echo "Normalized profile updated at: ${PROFILE}"

echo "-------------------------------------------------"
echo "Pipeline finished successfully for ${TARGET_RNA}."
echo "================================================="