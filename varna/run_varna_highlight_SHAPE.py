#!/usr/bin/env python3

# ==============================================================================
# Description:
#   Generates 2D RNA structure diagrams using VARNA. This script takes a summary
#   of RNA folding predictions and corresponding SHAPE reactivity data to produce
#   SVG images. Each nucleotide is colored based on its SHAPE reactivity (low,
#   medium, high). If a specific position of interest is provided for an RNA,
#   it is highlighted with a special color and a surrounding box. If no
#   highlight position is found, the structure is still generated without it.
#
# How to run:
#   ml devel python/3.9.0
#   source /path/to/your/python/venv/bin/activate
#   python3 /path/to/this/script.py --input-csv <path> --shape-dir <path> ...
#
# To create a log file:
#   python3 /path/to/this/script.py [ARGUMENTS] > varna_run.log 2>&1
# ==============================================================================

import pandas as pd
import os
import subprocess
import time
import argparse

# # --- Configuration ---
VARNA_JAR_PATH = os.path.expanduser('/home/users/rodell/varna/VARNAv3-93.jar')
# # 1. Path to your main CSV file with structures, from fold_summary.R
# INPUT_FILE_PATH = '/scratch/groups/nicolemm/rodell/SHAPE_InVitro/Pool1_Illumina/NoPUS_down/NoPUS_down_rnafold_summary.csv'  # <--- UPDATE THIS TO YOUR FILENAME
# # 2. Path to the directory containing individual .shape files for each RNA.
# SHAPE_DIR = '/scratch/groups/nicolemm/rodell/SHAPE_InVitro/Pool1_Illumina/NoPUS_down/shape_files'
# # 3. Path to your CSV defining the highlight position for each RNA.
# #    Must contain 'name' and 'psipos_noadapt' columns, where psipos_noadapt will be highlighted.
# HIGHLIGHT_CSV_PATH = '/home/groups/nicolemm/rodell/pool1/pool1_cleaned_20251107.csv'
# # HIGHLIGHT_POSITION = 66
# OUTPUT_DIR = "illumina_varna_structures_svg"

# --- Argument Parsing ---
parser = argparse.ArgumentParser(
    description="Generate VARNA 2D structure diagrams with SHAPE-based coloring and an optional highlighted position."
)
parser.add_argument("--input-csv", required=True, help="Path to the main CSV file containing sequence and structure information (e.g., from RNAfold).")
parser.add_argument("--shape-dir", required=True, help="Path to the directory containing per-RNA .shape files.")
parser.add_argument("--highlight-csv", required=True, help="Path to the CSV file defining the position to highlight for each RNA. Must contain columns for the RNA name and the position.")
parser.add_argument("--output-dir", required=True, help="Path to the directory where output SVG files will be saved.")
#parser.add_argument("--varna-jar", default="os.path.expanduser('/home/users/rodell/varna/VARNAv3-93.jar')", help="Path to the VARNAv3-93.jar executable file. Default: os.path.expanduser('/home/users/rodell/varna/VARNAv3-93.jar')")
parser.add_argument("--rna-name-col", default="sequence_id", help="Name of the column in the input CSV that contains the RNA identifier. Default: 'sequence_id'.")
parser.add_argument("--highlight-name-col", default="clean_name", help="Name of the column in the highlight CSV that contains the RNA identifier. Default: 'clean_name'.")
parser.add_argument("--highlight-pos-col", default="psipos_noadapt", help="Name of the column in the highlight CSV that contains the position to highlight. Default: 'psipos_noadapt'.")

args = parser.parse_args()

# --- Print Configuration and Create Output Directory ---
print("--- VARNA High-Throughput Structure Generation ---")
print(f"Input Structure CSV: {args.input_csv}")
print(f"SHAPE Directory:     {args.shape_dir}")
print(f"Highlight CSV:       {args.highlight_csv}")
print(f"Output Directory:    {args.output_dir}")
#print(f"VARNA Jar Path:      {args.varna_jar}")
print("-" * 50)

os.makedirs(args.output_dir, exist_ok=True)

print("--- VARNA High-Throughput Structure Generation ---")

# 1. Load the DataFrame
try:
    df = pd.read_csv(args.input_csv) 
    highlight_df = pd.read_csv(args.highlight_csv)
    print(f"Successfully loaded {len(df)} records from input CSV.")
    print(f"Successfully loaded {len(highlight_df)} highlight positions from highlight CSV.")
except FileNotFoundError:
    print(f"Error: A required file was not found. Please check paths. Details: {e}")
    exit()

# Merge the highlight positions into the main dataframe
df = pd.merge(df, highlight_df, left_on=args.rna_name_col, right_on=args.highlight_name_col, how='left')
df[args.highlight_pos_col] = df[args.highlight_pos_col].astype('Int64')

df['low_reactivity_bases'] = ''
df['medium_reactivity_bases'] = ''
df['high_reactivity_bases'] = ''

# Iterate through the DataFrame to add SHAPE data
for index, row in df.iterrows():
    rna_name = row[args.rna_name_col]
    shape_filepath = os.path.join(args.shape_dir, f"{rna_name}.shape")

    if not os.path.exists(shape_filepath):
        print(f"  - Warning: No .shape file found for {rna_name}. Skipping SHAPE coloring.")
        continue

    # Read shape data and bin reactivities
    shape_data = pd.read_csv(shape_filepath, sep='\s+', header=None, names=['position', 'reactivity'])
    
    low_bases, medium_bases, high_bases = [], [], []
    for _, shape_row in shape_data.iterrows():
        base_num_str = str(int(shape_row['position']))
        reactivity = shape_row['reactivity']
        
        if reactivity < 0.3:
            low_bases.append(base_num_str)
        elif reactivity <= 0.7:
            medium_bases.append(base_num_str)
        else:
            high_bases.append(base_num_str)
            
    # Store the binned lists (as comma-separated strings for VARNA)
    df.at[index, 'low_reactivity_bases'] = ",".join(low_bases)
    df.at[index, 'medium_reactivity_bases'] = ",".join(medium_bases)
    df.at[index, 'high_reactivity_bases'] = ",".join(high_bases)

print("--- Data preparation complete. ---")

# --- Main Loop: Iterate over every RNA in the DataFrame ---
total_rnas = len(df)
start_time = time.time()

for index, row in df.iterrows():
    rna_name = row[args.rna_name_col]
    highlight_pos = row[args.highlight_pos_col]
    print(f"\nProcessing RNA {index + 1}/{total_rnas}: {rna_name}")

    # Prepare file paths
    safe_filename = rna_name.replace(',', '_').replace('/', '_')
    structure_filepath = os.path.join(args.output_dir, f"{safe_filename}.db")
    output_filepath = os.path.join(args.output_dir, f"{safe_filename}.svg")

    # Create the temporary .db file for this RNA
    with open(structure_filepath, "w") as f:
        f.write(f">{rna_name}\n")
        f.write(f"{row['sequence']}\n") 
        f.write(f"{row['mea_structure']}\n")

    # --- Construct the VARNA command with your final theme ---
    command = [
        "java", "-cp", VARNA_JAR_PATH, "fr.orsay.lri.varna.applications.VARNAcmd",
        "-i", structure_filepath,
        "-o", output_filepath,
        # General options
        "-title", rna_name, "-resolution", "2.0", "-periodNum", "10",
        # Style 1: Low reactivity (grey)
        "-basesStyle1", "outline=#000000,label=#000000,fill=#e5e4e2",
        # Style 2: Medium reactivity (orange)
        "-basesStyle2", "outline=#b26200,label=#b26200,fill=#ffc87c",
        # Style 3: High reactivity (red)
        "-basesStyle3", "outline=#6b1010,label=#6b1010,fill=#e0a3a3",
        # Style 4: Special highlight position (pink label)
        "-basesStyle4", "label=#c63e83,outline=#000000",
    ]

    # --- Apply styles based on data ---
    # Apply general SHAPE styles first
    if row['low_reactivity_bases']:
        command.extend(["-applyBasesStyle1on", row['low_reactivity_bases']])
    if row['medium_reactivity_bases']:
        command.extend(["-applyBasesStyle2on", row['medium_reactivity_bases']])
    if row['high_reactivity_bases']:
        command.extend(["-applyBasesStyle3on", row['high_reactivity_bases']])

    if pd.notna(highlight_pos):
        print("  - Highlight position found. Applying special styles.")
        command.extend(["-applyBasesStyle4on", str(highlight_pos)])
        highlight_region_value = f"{highlight_pos}-{highlight_pos}:fill=#edc5d9,outline=#c63e83,radius=25"
        command.extend(["-highlightRegion", highlight_region_value])
    else:
        print("  - No highlight position found. Proceeding without special highlight.")

    # --- Execute the command for this RNA ---
    # print(command)
    try:
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        print(f"  - Successfully created image: {output_filepath}")
    except subprocess.CalledProcessError as e:
        print(f"  - !!! ERROR on {rna_name} !!!")
        print("    VARNA command failed to execute. Skipping this RNA.")
        print(f"    Return Code: {e.returncode}\n    VARNA Error Output (stderr):\n    {e.stderr}")
    
    # Clean up the temporary .db file to save space
    os.remove(structure_filepath)

# --- Final Summary ---
end_time = time.time()
duration = end_time - start_time
print("\n-------------------------------------------------")
print("High-throughput processing complete.")
print(f"Processed {total_rnas} RNAs in {duration:.2f} seconds.")
print(f"Output saved to '{args.output_dir}'.")
print("-------------------------------------------------")
