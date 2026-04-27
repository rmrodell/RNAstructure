#!/bin/bash

ml biology viennarna

echo "Input .fa file: $1"
echo "Output directory: $2"
echo "Current working directory: $(pwd)"

input_file=$(realpath "$1")
echo "Input file real path: $input_file"

# Check if input file exists
if [ ! -f "$1" ]; then
    echo "Error: Input file $1 does not exist"
    exit 1
fi

# Check if output directory exists, if not create it
if [ ! -d "$2" ]; then
    mkdir -p "$2"
fi

# Get the base name of the input file
base_name=$(basename "$1" .fa)

# Change to the output directory
cd "$2"

# Run RNAfold
echo "Running RNAfold..."
RNAfold -p -o --MEA "$input_file"


# Check if RNAfold was successful
if [ $? -eq 0 ]; then
    echo "RNAfold completed successfully"
else
    echo "Error: RNAfold failed"
    exit 1
fi

echo "RNAfold output files saved in $2"