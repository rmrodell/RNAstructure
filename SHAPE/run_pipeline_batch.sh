#!/bin/bash
#
# Wrapper script to run the SHAPE-MaP pipeline for multiple samples in batch.
# It reads a sample map file and launches a separate pipeline run for each entry.
#
# Usage: ./batch_shape_pipeline.sh \
#         --bam_source_dir <bam_source_dir> \
#         --sample_map_file <sample_map_file> \
#         --top_level_output_dir <top_level_output_dir> \
#         --pipeline_script <pipeline_script>
#
# Example: ./batch_shape_pipeline.sh \
#         --bam_source_dir /scratch/groups/nicolemm/rodell/basecalling/InVitro_SHAPE_Rep1_20241209/final_results/ \
#         --sample_map_file /scratch/groups/nicolemm/rodell/SHAPE/InVitro_Rep1_20241209/samples.txt \
#         --top_level_output_dir /scratch/groups/nicolemm/rodell/SHAPE/InVitro_Rep1_20241209/batch_output_test \
#         --pipeline_script /scratch/groups/nicolemm/rodell/SHAPE/process_single_file.sh


# Stop the script if any command fails
set -e
set -o pipefail

################################################################################
### USER-CONFIGURABLE VARIABLES                                            ###
################################################################################

# Initialize variables with default values (if desired)
BAM_SOURCE_DIR=""
SAMPLE_MAP_FILE=""
TOP_LEVEL_OUTPUT_DIR=""
PIPELINE_SCRIPT="/scratch/groups/nicolemm/rodell/SHAPE/process_single_file.sh"

# Function to display usage information
usage() {
  echo "Usage: $0 --bam_source_dir <bam_source_dir> --sample_map_file <sample_map_file> --top_level_output_dir <top_level_output_dir> --pipeline_script <pipeline_script>"
  echo "  --bam_source_dir         Directory containing the raw input BAM files."
  echo "  --sample_map_file        File containing the 'SampleName:barcode' mapping."
  echo "  --top_level_output_dir   Parent directory for all outputs."
  echo "  --pipeline_script        Path to the pipeline script to execute for each sample."
  exit 1
}

# Parse command-line options
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam_source_dir)
      BAM_SOURCE_DIR="$2"
      shift 2
      ;;
    --sample_map_file)
      SAMPLE_MAP_FILE="$2"
      shift 2
      ;;
    --top_level_output_dir)
      TOP_LEVEL_OUTPUT_DIR="$2"
      shift 2
      ;;
    --pipeline_script)
      PIPELINE_SCRIPT="$2"
      shift 2
      ;;
    *)
      echo "Unknown parameter: $1"
      usage
      ;;
  esac
done

# Check if all required options are provided
if [ -z "$BAM_SOURCE_DIR" ] || [ -z "$SAMPLE_MAP_FILE" ] || [ -z "$TOP_LEVEL_OUTPUT_DIR" ] || [ -z "$PIPELINE_SCRIPT" ]; then
  echo "Error: Missing required parameters."
  usage
fi

echo "BAM_SOURCE_DIR: $BAM_SOURCE_DIR"
echo "SAMPLE_MAP_FILE: $SAMPLE_MAP_FILE"
echo "TOP_LEVEL_OUTPUT_DIR: $TOP_LEVEL_OUTPUT_DIR"
echo "PIPELINE_SCRIPT: $PIPELINE_SCRIPT"

################################################################################
### BATCH EXECUTION LOGIC                                                  ###
################################################################################

echo "--- Starting Batch Pipeline Run ---"

# Define the common directory where final FASTQ files will be placed.
FINAL_FASTQ_DIR="${TOP_LEVEL_OUTPUT_DIR}/final_deduplicated_fastq"

# Ensure the main output directory and the final results subdirectory exist
mkdir -p "$TOP_LEVEL_OUTPUT_DIR"
mkdir -p "FINAL_FASTQ_DIR"
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
    "$PIPELINE_SCRIPT" "$BAM_FILE_PATH" "$sample_name" "$TOP_LEVEL_OUTPUT_DIR"

    echo "  - Finished processing for sample: $sample_name"

done < "$SAMPLE_MAP_FILE"

echo "-------------------------------------------------------------"
echo "--- Batch run completed successfully! ---"