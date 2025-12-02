#!/bin/bash
#
# 01_shapemapper_all_targets.sh
#
# Wrapper script to iterate over RNA targets from a FASTA file and run the
# shapemapper_single_target.sh pipeline for each.

set -e
set -o pipefail

# --- Usage function ---
usage() {
    echo "Usage: $0 -d <shapemapper_dir> -f <fasta_file> -m <modified_bam> -u <untreated_bam> -o <output_dir> -p <prefix>"
    echo "  -d  Path to the main shapemapper2-2.3 directory."
    echo "  -f  Path to the multi-sequence FASTA reference file."
    echo "  -m  Path to the 'modified' sample sorted BAM file."
    echo "  -u  Path to the 'untreated' (e.g., DMSO) sample sorted BAM file."
    echo "  -o  Path to the main output directory for all results."
    echo "  -p  A prefix for output file names (e.g., 'Pool2_Rep1')."
    exit 1
}

# --- Parse Command-Line Arguments ---
# Note: No -r argument, as we will get it from the FASTA file
while getopts ":d:f:m:u:o:p:" opt; do
    case ${opt} in
        d) SHAPEMAPPER_DIR=$OPTARG ;;
        f) FASTA_FILE=$OPTARG ;;
        m) MODIFIED_BAM=$OPTARG ;;
        u) UNTREATED_BAM=$OPTARG ;;
        o) OUTPUT_DIR=$OPTARG ;;
        p) PREFIX=$OPTARG ;;
        \?) echo "Invalid option: -$OPTARG" >&2; usage ;;
        :) echo "Option -$OPTARG requires an argument." >&2; usage ;;
    esac
done

# --- Check for completeness of arguments ---
if [ -z "$SHAPEMAPPER_DIR" ] || [ -z "$FASTA_FILE" ] || [ -z "$MODIFIED_BAM" ] || [ -z "$UNTREATED_BAM" ] || [ -z "$OUTPUT_DIR" ] || [ -z "$PREFIX" ]; then
    echo "Error: Missing one or more required arguments."
    usage
fi

# --- Find the single-target script ---
SCRIPT_DIR=$(dirname "$0")
SINGLE_TARGET_SCRIPT="${SCRIPT_DIR}/shapemapper_single_target.sh"

if [ ! -f "$SINGLE_TARGET_SCRIPT" ]; then
    echo "Error: The single target script was not found at ${SINGLE_TARGET_SCRIPT}"
    exit 1
fi

# # --- Setup Directories and Logging for THIS script ---
mkdir -p "$OUTPUT_DIR"
LOG_FILE="${OUTPUT_DIR}/wrapper_pipeline_run.log"
exec > "$LOG_FILE" 2>&1

# --- Define path for the list of failed targets and initialize it ---
FAILED_TARGETS_LIST="${OUTPUT_DIR}/${PREFIX}_failed_targets.txt"
> "$FAILED_TARGETS_LIST" # Create or overwrite the file to ensure it's empty

# --- Main Loop ---
echo "================================================="
echo "Starting ShapeMapper Pipeline for All Targets"
echo "================================================="

# Extract RNA names from the FASTA file.
RNA_TARGETS=$(grep '^>' "$FASTA_FILE" | awk '{print substr($1,2)}')

# For testing with just 15 sequences:
# RNA_TARGETS=$(echo "$RNA_TARGETS" | head -n 15)


echo "Found the following targets to process (first 10 shown):"
echo "$RNA_TARGETS" | head -n 10
echo "-------------------------------------------------"

# Loop over each target and call the single-target script
for target in $RNA_TARGETS; do
    echo ">>>> Processing target: ${target} <<<<"

    if "$SINGLE_TARGET_SCRIPT" \
        -d "$SHAPEMAPPER_DIR" \
        -f "$FASTA_FILE" \
        -r "$target" \
        -m "$MODIFIED_BAM" \
        -u "$UNTREATED_BAM" \
        -o "$OUTPUT_DIR" \
        -p "$PREFIX"
    then
        echo "--- Finished processing ${target}. Cleaning up intermediate files. ---"
        
        # Define intermediate files to delete
        MOD_SAM_FILTERED="${OUTPUT_DIR}/${PREFIX}_modified_${target}.sam"
        UN_SAM_FILTERED="${OUTPUT_DIR}/${PREFIX}_untreated_${target}.sam"
        MOD_MUT="${OUTPUT_DIR}/${PREFIX}_modified_${target}.mut"
        UN_MUT="${OUTPUT_DIR}/${PREFIX}_untreated_${target}.mut"
        MOD_COUNTS="${OUTPUT_DIR}/${PREFIX}_modified_${target}_counts.txt"
        UN_COUNTS="${OUTPUT_DIR}/${PREFIX}_untreated_${target}_counts.txt"
        SINGLE_TARGET_LOG="${OUTPUT_DIR}/${PREFIX}_${target}_pipeline.log"
        
        # Delete the files, leaving only the final *_profile.txt
        rm -f \
            "$MOD_SAM_FILTERED" \
            "$UN_SAM_FILTERED" \
            "$MOD_MUT" \
            "$UN_MUT" \
            "$MOD_COUNTS" \
            "$UN_COUNTS" \
            "$SINGLE_TARGET_LOG"
            
        echo "--- Cleanup complete. Final profile retained: ${OUTPUT_DIR}/${PREFIX}_${target}_profile.txt ---"
        echo ""
    else
        # This block executes if the single_target_script fails (returns a non-zero exit code)
        echo "!!!!!! FAILURE: Target ${target} failed. Appending to list. !!!!!!"
        echo "$target" >> "$FAILED_TARGETS_LIST"
        echo "!!!!!! Intermediate files are being kept for debugging. !!!!!! "
    fi
    
    echo ""
done

echo "================================================="
echo "All targets have been processed."
echo "================================================="