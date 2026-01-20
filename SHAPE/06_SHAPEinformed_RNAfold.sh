#!/bin/bash

# =========================================================================
# Description:
#   This script automates running ViennaRNA's RNAfold with SHAPE reactivity data.
#   It iterates through each sequence in a multi-FASTA file, finds its corresponding
#   .shape file in a specified directory, and executes RNAfold.
#
#   All output from the script (progress, warnings, errors) is saved to a log file.
#
# Arguments:
#   -f, --fasta      <path>   Path to the reference FASTA file. (Required)
#   -s, --shape-dir  <path>   Path to the directory containing the .shape files.
#                             (Required)
#   -o, --out-dir    <path>   Path to the directory where all RNAfold output files
#                             will be saved. (Required)
#   -h, --help                Display this help message and exit.
#
# Module Requirements:
#   This script requires 'viennarna' to be available. On an HPC, this is typically
#   handled by a module system (e.g., 'module load viennarna').
#
# Example Usage:
#   ./SHAPEinformed_RNAfold.sh \
#       --fasta /path/to/sequences.fasta \
#       --shape-dir /path/to/shape_files \
#       --out-dir /path/to/rnafold_output
# =========================================================================

# --- Define a function for displaying usage information ---
usage() {
    echo "Usage: $0 -f <fasta_file> -s <shape_dir> -o <out_dir>"
    echo ""
    echo "Required Arguments:"
    echo "  -f, --fasta      <path>   Path to the reference FASTA file."
    echo "  -s, --shape-dir  <path>   Directory containing corresponding .shape files."
    echo "  -o, --out-dir    <path>   Directory to save RNAfold outputs."
    echo ""
    echo "Optional Arguments:"
    echo "  -h, --help                Display this help message and exit."
    exit 1
}

# --- Initialize variables ---
FASTA_FILE=""
SHAPE_DIR=""
OUTPUT_DIR=""

# --- 1. USER-DEFINED PATHS ---
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -f|--fasta) FASTA_FILE="$2"; shift ;;
        -s|--shape-dir) SHAPE_DIR="$2"; shift ;;
        -o|--out-dir) OUTPUT_DIR="$2"; shift ;;
        -h|--help) usage ;;
        *) echo "ERROR: Unknown parameter passed: $1"; usage ;;
    esac
    shift
done

# Check for missing required arguments
if [[ -z "$FASTA_FILE" || -z "$SHAPE_DIR" || -z "$OUTPUT_DIR" ]]; then
    echo "ERROR: One or more required arguments are missing."
    usage
fi

# Validate that inputs exist
if [[ ! -f "$FASTA_FILE" ]]; then
    echo "ERROR: FASTA file not found at: ${FASTA_FILE}"
    exit 1
fi
if [[ ! -d "$SHAPE_DIR" ]]; then
    echo "ERROR: SHAPE directory not found at: ${SHAPE_DIR}"
    exit 1
fi

# --- 2. SETUP AND LOGGING ---
mkdir -p "$OUTPUT_DIR"
LOG_DIR="${OUTPUT_DIR}/logs"
mkdir -p "$LOG_DIR"
LOG_FILE="${LOG_DIR}/rnafold_run_$(date +%Y-%m-%d_%H-%M-%S).log"

# Redirect all subsequent script output (stdout and stderr) only to the log file.
exec &> "$LOG_FILE"

echo "===================================================="
echo "Starting SHAPE-informed RNAfold analysis on $(date)"
echo "----------------------------------------------------"
echo "Parameters provided:"
echo "  - FASTA File:      ${FASTA_FILE}"
echo "  - SHAPE Directory: ${SHAPE_DIR}"
echo "  - Output Directory:  ${OUTPUT_DIR}"
echo ""
echo "Configuration:"
echo "  - Log File:        ${LOG_FILE}"
echo "===================================================="

# Load necessary software modules for the cluster environment.
ml biology viennarna

# --- 3. EXECUTION LOGIC ---

# Create a clean list of RNA names from the FASTA file headers.
RNA_NAME_LIST=$(grep '^>' "$FASTA_FILE" | sed 's/>//' | awk '{print $1}')

# Change to the output directory. All subsequent commands will run from here,
# and all output files will be written here directly.
cd "$OUTPUT_DIR"

# Loop through each clean name in the list.
for rna_name in $RNA_NAME_LIST; do
    
    # Sanitize the name for filename matching by replacing commas with underscores.
    sanitized_rna_name=$(echo "$rna_name" | tr ',' '_')
    shape_file="${SHAPE_DIR}/${sanitized_rna_name}.shape"

    # Check if the corresponding .shape file exists.
    if [ -f "$shape_file" ]; then
        echo "Processing: $rna_name"

        # Extract the sequence for this specific RNA.
        # This command finds the header line, gets the next line, and joins it.
        sequence=$(awk -v name="^>${rna_name}$" 'BEGIN{p=0} $0~name{p=1;next} /^>/{p=0} p{printf "%s", $0}' "$FASTA_FILE")

        # Run RNAfold. Output files (e.g., RNA_NAME_ss.ps) will be created in the current directory.
        RNAfold -p --shape="$shape_file" --MEA <<< ">${rna_name}
${sequence}" > "${rna_name}.fold"

    else
        echo "WARNING: Skipping '$rna_name'. No matching .shape file found."
    fi
done

echo "----------------------------------------------------"
echo "RNAfold analysis complete. All output is in ${OUTPUT_DIR}"