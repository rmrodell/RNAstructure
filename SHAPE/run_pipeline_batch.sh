#!/bin/bash
#
# Wrapper script to run the SHAPE-MaP pipeline for multiple samples in batch.
# It reads a sample map file and launches a separate pipeline run for each entry.
#

# Stop the script if any command fails
set -e
set -o pipefail

################################################################################
### USER-CONFIGURABLE VARIABLES                                            ###
################################################################################

# Directory containing the raw input BAM files from basecalling.
BAM_SOURCE_DIR="/scratch/groups/nicolemm/rodell/basecalling/InVitro_SHAPE_Rep1_20241209/final_results/"

# The file containing the "SampleName:barcode" mapping.
SAMPLE_MAP_FILE="/scratch/groups/nicolemm/rodell/SHAPE/InVitro_Rep1_20241209/samples.txt"

# --- The SINGLE parent directory for ALL outputs (intermediate files and final results) ---
TOP_LEVEL_OUTPUT_DIR="/scratch/groups/nicolemm/rodell/SHAPE/InVitro_Rep1_20241209/batch_output_test"

# The location of the pipeline script that will be called for each sample.
PIPELINE_SCRIPT="/scratch/groups/nicolemm/rodell/SHAPE/process_single_file.sh"

################################################################################
### BATCH EXECUTION LOGIC                                                  ###
################################################################################

echo "--- Starting Batch Pipeline Run ---"

# Define the common directory where final FASTQ files will be placed.
FINAL_FASTQ_DIR="${TOP_LEVEL_OUTPUT_DIR}/final_deduplicated_fastq"

# Ensure the main output directory and the final results subdirectory exist
mkdir -p "$TOP_LEVEL_OUTPUT_DIR"
mkdir -p "$FINAL_OUTPUT_DIR"
echo "All outputs will be placed within: $TOP_LEVEL_OUTPUT_DIR"
echo "Final FASTQ files will be copied to: $FINAL_FASTQ_DIR"

# Check if the pipeline script is executable
if [ ! -x "$PIPELINE_SCRIPT" ]; then
    echo "Error: Pipeline script '$PIPELINE_SCRIPT' not found or not executable."
    echo "Please run: chmod +x $PIPELINE_SCRIPT"
    exit 1
fi

# Read the sample map file line by line
while IFS=':' read -r sample_name barcode || [[ -n "$barcode" ]]; do
    # Trim whitespace just in case
    sample_name=$(echo "$sample_name" | xargs)
    barcode=$(echo "$barcode" | xargs)

    if [ -z "$sample_name" ] || [ -z "$barcode" ]; then
        echo "Skipping invalid line in $SAMPLE_MAP_FILE..."
        continue
    fi

    echo "-------------------------------------------------------------"
    echo "Processing Sample: $sample_name (Barcode: $barcode)"

    # Find the corresponding BAM file
    # This assumes the file format is *_{barcode}.bam
    BAM_FILE_PATH=$(find "$BAM_SOURCE_DIR" -name "*_${barcode}.bam")

    if [ -z "$BAM_FILE_PATH" ]; then
        echo "  - WARNING: BAM file for barcode '$barcode' not found in '$BAM_SOURCE_DIR'. Skipping."
        continue
    elif [ $(echo "$BAM_FILE_PATH" | wc -l) -gt 1 ]; then
        echo "  - WARNING: Multiple BAM files found for barcode '$barcode'. Skipping."
        continue
    fi

    echo "  - Found BAM file: $BAM_FILE_PATH"
    echo "  - Launching pipeline..."

    # Execute the main pipeline script, passing the required information as arguments
    "$PIPELINE_SCRIPT" "$BAM_FILE_PATH" "$sample_name" "$TOP_LEVEL_OUTPUT_DIR" "$FINAL_FASTQ_DIR"

    echo "  - Finished processing for sample: $sample_name"

done < "$SAMPLE_MAP_FILE"

echo "-------------------------------------------------------------"
echo "--- Batch run completed successfully! ---"