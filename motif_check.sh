#!/bin/bash

ml biology viennarna
ml R

out_dir_bestf1="/scratch/users/rodell/motifmatcher/20250802/analysis/bestf1"
out_dir_mostmut="/scratch/users/rodell/motifmatcher/20250802/analysis/mostmut"

mkdir -p $out_dir_bestf1
mkdir -p $out_dir_mostmut

# --- Run motif_matcher_v3.R with the given parameter sets. ---

# best-f1
Rscript /scratch/users/rodell/motifmatcher/RNAstructure/motif_matcher_v3.R \
    --input "/scratch/users/rodell/motifmatcher/input_pool1.csv" --fold_dir "/scratch/users/rodell/RNAfold_psipos" \
    --output "${out_dir_bestf1}/motifmatcher1.csv" \
    --offset_min 0 --offset_max 2 \
    --min_unpaired1 1 --max_unpaired1 4 \
    --min_paired1 2 --max_paired1 7 \
    --min_unpaired2 3 --max_unpaired2 5 \
    --min_paired2 1 --max_paired2 4 \
    --include_unpaired1 TRUE --include_paired2 TRUE

# most-mut
Rscript /scratch/users/rodell/motifmatcher/RNAstructure/motif_matcher_v3.R \
    --input "/scratch/users/rodell/motifmatcher/input_pool1.csv" --fold_dir "/scratch/users/rodell/RNAfold_psipos" \
    --output "${out_dir_mostmut}/motifmatcher1.csv" \
    --offset_min 0 --offset_max 2 \
    --min_unpaired1 1 --max_unpaired1 8 \
    --min_paired1 2 --max_paired1 9 \
    --min_unpaired2 3 --max_unpaired2 9 \
    --min_paired2 1 --max_paired2 4 \
    --include_unpaired1 TRUE --include_paired2 TRUE


# --- Perform mutagenesis with mutagenesis.R for the matching sequences. ---

# best-f1
Rscript /scratch/users/rodell/motifmatcher/RNAstructure/mutagenesis.R \
    --input "${out_dir_bestf1}/motifmatcher1.csv" --output "${out_dir_bestf1}/mutagenesis.csv" --protect_position 60

# most-mut
Rscript /scratch/users/rodell/motifmatcher/RNAstructure/mutagenesis.R \
    --input "${out_dir_mostmut}/motifmatcher1.csv" --output "${out_dir_mostmut}/mutagenesis.csv" --protect_position 60


# --- Write out mutated sequences to fasta with write_mutant_fastas.R ---

# best-f1
Rscript /scratch/users/rodell/motifmatcher/RNAstructure/write_mutant_fastas.R \
    "${out_dir_bestf1}/mutagenesis.csv" "${out_dir_bestf1}"

# most-mut
Rscript /scratch/users/rodell/motifmatcher/RNAstructure/write_mutant_fastas.R \
    "${out_dir_mostmut}/mutagenesis.csv" "${out_dir_mostmut}"



# --- Run RNAfold on mutated sequences with run_RNAfold.sh ---

SCRIPT="/scratch/users/rodell/motifmatcher/RNAstructure/run_RNAfold.sh"

# List of mutants to process
mutant_types=(
    "paired1_compensatory"
    "paired1_disruption"
    "paired2_compensatory"
    "paired2_disruption"
)

# best-f1
for mutant in "${mutant_types[@]}"; do
    output_dir="${out_dir_bestf1}/${mutant}"
    fasta_file="${output_dir}.fa"

    mkdir -p "$output_dir"
    bash "$SCRIPT" "$fasta_file" "$output_dir"
done

# most-mut
for mutant in "${mutant_types[@]}"; do
    output_dir="${out_dir_mostmut}/${mutant}"
    fasta_file="${output_dir}.fa"

    mkdir -p "$output_dir"
    bash "$SCRIPT" "$fasta_file" "$output_dir"
done

# Run motif_matcher_v3.R on mutated sequences.

# best-f1
for mutant in "${mutant_types[@]}"; do

    Rscript /scratch/users/rodell/motifmatcher/RNAstructure/motif_matcher_v3.R \
        --input "/scratch/users/rodell/motifmatcher/input_pool1.csv" --fold_dir "${out_dir_bestf1}/${mutant}" \
        --output "${out_dir_bestf1}/${mutant}_motifmatcher2.csv" \
        --offset_min 0 --offset_max 2 \
        --min_unpaired1 1 --max_unpaired1 4 \
        --min_paired1 2 --max_paired1 7 \
        --min_unpaired2 3 --max_unpaired2 5 \
        --min_paired2 1 --max_paired2 4 \
        --include_unpaired1 TRUE --include_paired2 TRUE

done

# most-mut
for mutant in "${mutant_types[@]}"; do

    Rscript /scratch/users/rodell/motifmatcher/RNAstructure/motif_matcher_v3.R \
        --input "/scratch/users/rodell/motifmatcher/input_pool1.csv" --fold_dir "${out_dir_mostmut}/${mutant}" \
        --output "${out_dir_mostmut}/${mutant}_motifmatcher2.csv" \
        --offset_min 0 --offset_max 2 \
        --min_unpaired1 1 --max_unpaired1 8 \
        --min_paired1 2 --max_paired1 9 \
        --min_unpaired2 3 --max_unpaired2 9 \
        --min_paired2 1 --max_paired2 4 \
        --include_unpaired1 TRUE --include_paired2 TRUE

done

# Tabulate how many were disrupted (no longer match motif) and then compensated (match motif, but with different sequence)
# download csv files, work in local RStudio