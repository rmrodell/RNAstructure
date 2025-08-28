#!/bin/bash

ml biology viennarna
ml R


SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)

# --- Define arguments ---
INPUT_FILE="/scratch/users/rodell/motifmatcher/input_pool1.csv"
FOLD_DIR="/scratch/users/rodell/20250821_poollength/all_extend_trimmed"
OUTPUT_DIR="bestf1"

mkdir -p "$OUTPUT_DIR"

INPUT_POS=66
OFFSET_MIN=0
OFFSET_MAX=1
MIN_UNPAIRED1=1
MAX_UNPAIRED1=4
MIN_PAIRED1=3
MAX_PAIRED1=5
MIN_UNPAIRED2=2
MAX_UNPAIRED2=5
MIN_PAIRED2=1
MAX_PAIRED2=4
INCLUDE_UNPAIRED1=TRUE
INCLUDE_PAIRED2=TRUE

# --- Run motif_matcher_v3.R with the given parameter sets. ---

# best-f1
Rscript "${SCRIPT_DIR}/parameter_sweep_MFE/motif_matcher_v3.R" \
    --input "$INPUT_FILE" --fold_dir "$FOLD_DIR" \
    --output "$OUTPUT_DIR/motifmatcher1.csv" \
    --input_position "$INPUT_POS" \
    --offset_min "$OFFSET_MIN" --offset_max "$OFFSET_MAX" \
    --min_unpaired1 "$MIN_UNPAIRED1" --max_unpaired1 "$MAX_UNPAIRED1" \
    --min_paired1 "$MIN_PAIRED1" --max_paired1 "$MAX_PAIRED1" \
    --min_unpaired2 "$MIN_UNPAIRED2" --max_unpaired2 "$MAX_UNPAIRED2" \
    --min_paired2 "$MIN_PAIRED2" --max_paired2 "$MAX_PAIRED2" \
    --include_unpaired1 "$INCLUDE_UNPAIRED1" --include_paired2 "$INCLUDE_PAIRED2"

# --- Perform mutagenesis with mutagenesis.R for the matching sequences. ---

# best-f1
Rscript ${SCRIPT_DIR}/parameter_sweep_MFE/mutagenesis_v2.R \
    --input "$OUTPUT_DIR/motifmatcher1.csv" --output "$OUTPUT_DIR/mutagenesis.csv" --protect_position 66



# --- Write out mutated sequences to fasta with write_mutant_fastas.R ---

# best-f1
Rscript ${SCRIPT_DIR}/parameter_sweep_MFE/write_mutant_fastas.R \
    "$OUTPUT_DIR/mutagenesis.csv" "$OUTPUT_DIR"

# --- Run RNAfold on mutated sequences with run_RNAfold.sh ---

SCRIPT="${SCRIPT_DIR}/parameter_sweep_MFE/run_RNAfold.sh"

# List of mutants to process
mutant_types=(
    "paired1_compensatory"
    "paired1_disruption"
    "paired2_compensatory"
    "paired2_disruption"
    "combined_compensatory"
    "combined_disruption"
)

# best-f1
for mutant in "${mutant_types[@]}"; do
    output_dir="$OUTPUT_DIR/${mutant}"
    fasta_file="$OUTPUT_DIR/${mutant}.fa"

    # Check if the fasta file actually exists before processing
    if [[ ! -f "$fasta_file" ]]; then
        echo "Error: Fasta file not found for mutant type '$mutant' at: $fasta_file"
        continue # Skip to the next mutant type
    fi

    mkdir -p "$OUTPUT_DIR"
    bash "$SCRIPT" "$fasta_file" "$output_dir"
done

# Run motif_matcher_v3.R on mutated sequences.

# best-f1
for mutant in "${mutant_types[@]}"; do

    Rscript "${SCRIPT_DIR}/parameter_sweep_MFE/motif_matcher_v3.R" \
        --input "$INPUT_FILE" --fold_dir "$OUTPUT_DIR/${mutant}" \
        --output "$OUTPUT_DIR/${mutant}_motifmatcher2.csv" \
    --input_position "$INPUT_POS" \
    --offset_min "$OFFSET_MIN" --offset_max "$OFFSET_MAX" \
    --min_unpaired1 "$MIN_UNPAIRED1" --max_unpaired1 "$MAX_UNPAIRED1" \
    --min_paired1 "$MIN_PAIRED1" --max_paired1 "$MAX_PAIRED1" \
    --min_unpaired2 "$MIN_UNPAIRED2" --max_unpaired2 "$MAX_UNPAIRED2" \
    --min_paired2 "$MIN_PAIRED2" --max_paired2 "$MAX_PAIRED2" \
    --include_unpaired1 "$INCLUDE_UNPAIRED1" --include_paired2 "$INCLUDE_PAIRED2"

done

# Tabulate how many were disrupted (no longer match motif) and then compensated (match motif, but with different sequence)

