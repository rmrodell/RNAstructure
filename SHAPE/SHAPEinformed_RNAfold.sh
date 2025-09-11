#!/bin/bash

# =========================================================================
#       RUN SHAPE-Informed RNAfold for a Pool of Sequences
# =========================================================================
# This script reads a multi-FASTA file, and for each sequence, finds its
# corresponding averaged .shape file and runs RNAfold.

#set -e # Exit immediately if a command exits with a non-zero status.

# Load necessary software modules for the cluster environment.
# Adjust these if your cluster uses different module names.
ml biology viennarna

# --- 1. USER-DEFINED PATHS ---
FASTA_FILE="/home/groups/nicolemm/rodell/pool1/pool1_prettyplease_noadapters.fasta"
SHAPE_DIR="/scratch/groups/nicolemm/rodell/SHAPE/InVitro_2rep_Analysis_fullpool/averaged_shape_files"
OUTPUT_DIR="/scratch/groups/nicolemm/rodell/SHAPE/InVitro_2rep_Analysis_fullpool/RNAfold_output"
LOG_DIR="/scratch/groups/nicolemm/rodell/SHAPE/InVitro_2rep_Analysis_fullpool/run_logs"


# --- 2. SETUP AND LOGGING ---
mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"
LOG_FILE="${LOG_DIR}/rnafold_run_$(date +%Y-%m-%d_%H-%M-%S).log"

# Redirect all subsequent script output (stdout and stderr) to the log file.
exec > "$LOG_FILE" 2>&1

echo "Starting SHAPE-informed RNAfold analysis."
echo "Log file for this run: $LOG_FILE"
echo "FASTA file: $FASTA_FILE"
echo "SHAPE directory: $SHAPE_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "----------------------------------------------------"

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
        sequence=$(awk -v name="^>${rna_name}" 'BEGIN{p=0} $0~name{p=1;next} /^>/{p=0} p{printf "%s", $0}' "$FASTA_FILE")

        # Run RNAfold. Output files (e.g., RNA_NAME_ss.ps) will be created in the current directory.
        RNAfold -p --shape="$shape_file" --MEA <<< ">${rna_name}
${sequence}" > "${rna_name}.fold"

    else
        echo "WARNING: Skipping '$rna_name'. No matching .shape file found."
    fi
done

echo "----------------------------------------------------"
echo "RNAfold analysis complete. All output is in $OUTPUT_DIR"