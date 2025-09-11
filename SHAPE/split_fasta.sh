#!/bin/bash

# Split fasta into max of 200 chr / sequences for SHAPE-Mapper to work properly

# --- Configuration ---
# Define the path to your complete, original FASTA file
PARENT_FASTA="/home/groups/nicolemm/rodell/pool1/pool1_prettyplease_noadapters.fasta"

# Define the directory where you want the split FASTA chunks to be saved
FASTA_DIR="/home/groups/nicolemm/rodell/pool1/fasta_chunks"

# Define the prefix for the output chunk filenames
PREFIX="pool1"

# --- Execution ---
# Create the output directory if it doesn't already exist
echo "Creating output directory: $FASTA_DIR"
mkdir -p "$FASTA_DIR"

# Run the awk command, now with the custom prefix
echo "Splitting $PARENT_FASTA into chunks named '${PREFIX}_chunk_N.fasta'..."
awk -v out_dir="$FASTA_DIR" -v prefix="$PREFIX" '
    BEGIN {
        n = 1
        # Construct the filename using the prefix
        file = out_dir "/" prefix "_chunk_" n ".fasta"
    }
    /^>/ {
        if (count >= 100) {
            close(file)
            n++
            count = 0
            # Construct the filename for the next chunk
            file = out_dir "/" prefix "_chunk_" n ".fasta"
        }
        count++
    }
    {
        print > file
    }
' "$PARENT_FASTA"

echo "Done. FASTA chunks have been created in $FASTA_DIR"

# Verify the new filenames
echo "First 5 files created:"
ls -l "$FASTA_DIR"