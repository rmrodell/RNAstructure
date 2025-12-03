#!/bin/bash
#
# run_modular_pipeline.sh
#
# Master script to orchestrate the full SHAPE-MaP analysis pipeline for a single replicate.
#
# Pipeline Steps:
# 1. Runs 01_shapemapper_all_targets.sh to generate individual profiles for all targets.
#    This step now also generates a list of any targets that failed processing.
# 2. Runs the combined R script (02_combine_and_annotate.R) to:
#    a. Aggregate all individual profiles into a single file.
#    b. Annotate the combined file with QC flags based on the list of failed targets.
# 3. Runs the comparison R script (pool1_comp.R) to generate plots comparing the
#    newly generated profile data against a previous dataset.

set -e
set -o pipefail

# --- Usage function ---
usage() {
    echo "Usage: $0 -d <shapemapper_dir> -f <fasta_file> -m <modified_bam> -u <untreated_bam> -p <prefix> -c <comparison_file> -o <output_dir>"
    echo "  -d  Path to the main shapemapper2-2.3 directory."
    echo "  -f  Path to the multi-sequence FASTA reference file."
    echo "  -m  Path to the 'modified' sample sorted BAM file for this run."
    echo "  -u  Path to the 'untreated' (DMSO) sample sorted BAM file for this run."
    echo "  -p  A prefix for output file names for this run (e.g., 'Pool2_Rep1')."
    echo "  -c  Path to the final annotated profile file to compare against (e.g., Pool1's output)."
    echo "  -o  Path to the main parent output directory for all results."
    exit 1
}

# --- Parse Command-Line Arguments ---
while getopts ":d:f:m:u:p:c:o:" opt; do
    case ${opt} in
        d) SHAPEMAPPER_DIR=$OPTARG ;;
        f) FASTA_FILE=$OPTARG ;;
        m) MODIFIED_BAM=$OPTARG ;;
        u) UNTREATED_BAM=$OPTARG ;;
        p) PREFIX=$OPTARG ;;
        c) COMPARISON_FILE=$OPTARG ;;
        o) PARENT_OUTPUT_DIR=$OPTARG ;;
        \?) echo "Invalid option: -$OPTARG" >&2; usage ;;
        :) echo "Option -$OPTARG requires an argument." >&2; usage ;;
    esac
done

# --- Check for completeness of arguments ---
if [ -z "$SHAPEMAPPER_DIR" ] || [ -z "$FASTA_FILE" ] || [ -z "$MODIFIED_BAM" ] || [ -z "$UNTREATED_BAM" ] || [ -z "$PREFIX" ] || [ -z "$COMPARISON_FILE" ] || [ -z "$PARENT_OUTPUT_DIR" ]; then
    echo "Error: Missing one or more required arguments."
    usage
fi

# --- Setup Directories and File Paths ---
SCRIPT_DIR=$(dirname "$0")
SHAPEMAPPER_OUTPUT_DIR="${PARENT_OUTPUT_DIR}/shapemapper_profiles"
COMBINED_OUTPUT_DIR="${PARENT_OUTPUT_DIR}/combined_analysis"
COMPARISON_OUTPUT_DIR="${PARENT_OUTPUT_DIR}/comparison_plots"

mkdir -p "$SHAPEMAPPER_OUTPUT_DIR" "$COMBINED_OUTPUT_DIR" "$COMPARISON_OUTPUT_DIR"

# Define final file paths
COMBINED_ANNOTATED_PROFILE="${COMBINED_OUTPUT_DIR}/${PREFIX}_profile_annotated.txt"
FAILED_TARGETS_LIST="${SHAPEMAPPER_OUTPUT_DIR}/${PREFIX}_failed_targets.txt"

# --- Main Pipeline ---
echo "================================================="
echo "Starting Full SHAPE-MaP Analysis Pipeline"
echo "Prefix for this run: ${PREFIX}"
echo "Date: $(date)"
echo "================================================="

# --- Step 1: Generate individual SHAPE profiles ---
echo "[Step 1] Generating individual profiles for all targets..."
bash "${SCRIPT_DIR}/01_shapemapper_all_targets.sh" \
    -d "$SHAPEMAPPER_DIR" \
    -f "$FASTA_FILE" \
    -m "$MODIFIED_BAM" \
    -u "$UNTREATED_BAM" \
    -p "$PREFIX" \
    -o "$SHAPEMAPPER_OUTPUT_DIR"
echo "--- Individual profile generation complete. ---"
echo ""

# --- Step 2: Combine profiles and annotate QC flags ---
echo "[Step 2] Combining profiles and annotating QC flags..."
ml R
Rscript "${SCRIPT_DIR}/02_combine_annotate.R" \
    --dir "$SHAPEMAPPER_OUTPUT_DIR" \
    --pattern "${PREFIX}_.*_profile\\.txt" \
    --prefix "$PREFIX" \
    --failed-list "$FAILED_TARGETS_LIST" \
    --outfile "$COMBINED_ANNOTATED_PROFILE"
echo "--- Combination and annotation complete. ---"
echo ""

# --- Step 3: Run comparison against previous data ---
echo "[Step 3] Comparing generated profiles against reference data..."
Rscript "${SCRIPT_DIR}/pool1_comp.R" \
    --pool1 "$COMPARISON_FILE" \
    --pool2 "$COMBINED_ANNOTATED_PROFILE" \
    --out-prefix "${COMPARISON_OUTPUT_DIR}/${PREFIX}_vs_Pool1"
echo "--- Comparison complete. ---"
echo ""

echo "================================================="
echo "Full pipeline finished successfully."
echo "Final annotated profile: ${COMBINED_ANNOTATED_PROFILE}"
echo "Final comparison plots in: ${COMPARISON_OUTPUT_DIR}"
echo "================================================="