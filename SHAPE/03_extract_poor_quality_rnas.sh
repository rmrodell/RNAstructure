#!/bin/bash

# ===================================================================================
# Description:
#   This script parses a single SHAPE-Mapper log file to find warnings about "possible
#   poor quality reactivity profiles". It extracts the names of the associated
#   RNAs, creates a unique list, and saves it to a specified output file.
#
# Arguments:
#   -l, --log-file <path>  Path to SHAPE-Mapper log.
#                          (Required)
#   -o, --outfile <path>   Path for the final output file containing the list
#                          of RNA names. (Required)
#   -h, --help             Display this help message and exit.
#
# Example Usage:
#   ./extract_poor_quality_rnas.sh \
#       --log-file /path/to/shapemapper/log \
#       --outfile /path/to/Rep1_poor_quality_rnas.txt
#
# ===================================================================================

# --- Define a function for displaying usage information ---
usage() {
    echo "Usage: $0 -l <LOG_FILE> -p <log_pattern> -o <outfile>"
    echo ""
    echo "Required Arguments:"
    echo "  -l, --log-file <path>   Path to SHAPE-Mapper log."
    echo "  -o, --outfile <path>    Path for the output list of RNA names."
    echo ""
    echo "Optional Arguments:"
    echo "  -h, --help              Display this help message and exit."
    exit 1
}

# --- Initialize variables ---
LOG_FILE=""
OUTFILE=""

# =========================================================================
# STEP 1: ARGUMENT PARSING AND VALIDATION
# =========================================================================
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -l|--log-file) LOG_FILE="$2"; shift ;;
        -o|--outfile) OUTFILE="$2"; shift ;;
        -h|--help) usage ;;
        *) echo "ERROR: Unknown parameter passed: $1"; usage ;;
    esac
    shift
done

# Check for missing required arguments
if [[ -z "$LOG_FILE" || -z "$OUTFILE" ]]; then
    echo "ERROR: One or more required arguments are missing."
    usage
fi

# Check if the log directory exists
if [[ ! -f "$LOG_FILE" ]]; then
    echo "ERROR: Log file not found at: ${LOG_FILE}"
    exit 1
fi

# =========================================================================
# STEP 2: PREPARE ENVIRONMENT
# =========================================================================

# Ensure the output directory exists before trying to write to it
OUTPUT_DIR=$(dirname "${OUTFILE}")
mkdir -p "${OUTPUT_DIR}"
echo "Output will be saved to: ${OUTFILE}"

# =========================================================================
# STEP 3: EXECUTION LOGIC
# =========================================================================

echo "Processing log file: ${LOG_FILE}"

# Execute the extraction pipeline
# 1. grep: Find all lines with the warning text in the single log file.
# 2. awk #1: Isolate the comma-separated list of RNA names at the end of the line.
# 3. awk #2: Convert the comma-separated list into a newline-separated list.
# 4. sort -u: Sort the list and remove duplicate entries.
# 5. >: Redirect the final, unique list to the output file.
grep "possible poor quality reactivity profiles" "${LOG_FILE}" | \
    awk -F':' '{sub(/^[ ]+/, "", $NF); print $NF}' | \
    awk '{gsub(/, /,"\n"); print}' | \
    sort -u > "${OUTFILE}"

# Report the final count, handling the case where no matches are found
if [ -s "${OUTFILE}" ]; then
    RNA_COUNT=$(wc -l < "${OUTFILE}")
else
    RNA_COUNT=0
fi

echo "----------------------------------------------------"
echo "Processing complete."
echo "Found ${RNA_COUNT} unique RNAs with possible poor quality profiles."
echo "List saved to: ${OUTFILE}"