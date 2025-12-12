import pandas as pd
import os
import subprocess
import time

# How to run code:
# ml devel python/3.9.0
# source $HOME/RNAstructure/varna/varna-env/bin/activate
# python3 $HOME/RNAstructure/varna/run_varna_highlight.py

# --- Configuration ---
VARNA_JAR_PATH = os.path.expanduser('/home/users/rodell/varna/VARNAv3-93.jar')
# 1. Path to your main CSV file with structures, from fold_summary.R
INPUT_FILE_PATH = '/scratch/groups/nicolemm/rodell/SHAPE_InVitro/Pool1_Illumina/NoPUS_down/NoPUS_down_rnafold_summary.csv'  # <--- UPDATE THIS TO YOUR FILENAME
# 2. Path to the directory containing individual .shape files for each RNA.
SHAPE_DIR = '/scratch/groups/nicolemm/rodell/SHAPE_InVitro/Pool1_Illumina/NoPUS_down/shape_files'
# 3. Path to your CSV defining the highlight position for each RNA.
#    Must contain 'name' and 'psipos_noadapt' columns, where psipos_noadapt will be highlighted.
HIGHLIGHT_CSV_PATH = '/home/groups/nicolemm/rodell/pool1/pool1_cleaned_20251107.csv'
# HIGHLIGHT_POSITION = 66
OUTPUT_DIR = "illumina_varna_structures_svg"
os.makedirs(OUTPUT_DIR, exist_ok=True)


print("--- VARNA High-Throughput Structure Generation ---")

# 1. Load the DataFrame
try:
    df = pd.read_csv(INPUT_FILE_PATH) 
    highlight_df = pd.read_csv(HIGHLIGHT_CSV_PATH)
    print(f"Successfully loaded CSV file from '{INPUT_FILE_PATH}'.")
    print(f"Successfully loaded {len(highlight_df)} highlight positions from '{HIGHLIGHT_CSV_PATH}'.")
except FileNotFoundError:
    print(f"Error: A required file was not found. Please check paths. Details: {e}")
    exit()

total_rnas = len(df)
start_time = time.time()

# Merge the highlight positions into the main dataframe
# 'name' in highlight_df must match 'sequence_id' in df
df = pd.merge(df, highlight_df, left_on='sequence_id', right_on='clean_name', how='left')

# Initialize new columns to store prepared data
df['low_reactivity_bases'] = None
df['medium_reactivity_bases'] = None
df['high_reactivity_bases'] = None

# Iterate through the DataFrame to add SHAPE data
for index, row in df.iterrows():
    rna_name = row['sequence_id']
    shape_filepath = os.path.join(SHAPE_DIR, f"{rna_name}.shape")

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

for index, row in df.iterrows():
    rna_name = row['sequence_id']
    highlight_pos = row['psipos_noadapt']
    print(f"\nProcessing RNA {index + 1}/{total_rnas}: {rna_name}")

    # Skip if there is no highlight position defined for this RNA
    if pd.isna(highlight_pos):
        print("  - Warning: No highlight position found. Skipping this RNA.")
        continue

    # Prepare file paths
    safe_filename = rna_name.replace(',', '_').replace('/', '_')
    structure_filepath = os.path.join(OUTPUT_DIR, f"{safe_filename}.db")
    output_filepath = os.path.join(OUTPUT_DIR, f"{safe_filename}.svg")

    # Create the temporary .db file for this RNA
    with open(structure_filepath, "w") as f:
        f.write(f">{rna_name}\n")
        f.write(f"{row['sequence']}\n") 
        f.write(f"{row['mea_structure']}\n")

    # --- Construct the VARNA command with your final theme ---
    command = [
        # "java", "-cp", VARNA_JAR_PATH, "fr.orsay.lri.varna.applications.VARNAcmd",
        # "-i", structure_filepath, 
        # "-o", output_filepath,
        
        # # highlight residue style
        # "-basesStyle2", "fill=#dcdcdc,label=#c63e83,outline=#000000",
        # "-applyBasesStyle2on", str(HIGHLIGHT_POSITION),

        # # default style
        # "-basesStyle1", "outline=#000000,label=#000000,fill=#e5e4e2",
        # "-applyBasesStyle1on", "1-1000",

        # # General options
        # "-title", rna_name, "-resolution", "2.0", "-periodNum", "10",

        # "-highlightRegion", highlight_region_value

        "java", "-cp", VARNA_JAR_PATH, "fr.orsay.lri.varna.applications.VARNAcmd",
        "-i", structure_filepath,
        "-o", output_filepath,

        # General options
        "-title", rna_name, "-resolution", "2.0", "-periodNum", "10",

        # --- Define all styles ---
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
    
    # Apply the special highlight style LAST to ensure it overrides the SHAPE color
    command.extend(["-applyBasesStyle4on", str(int(highlight_pos))])

    # Add the highlight region box around the special position
    highlight_region_value = f"{int(highlight_pos)}-{int(highlight_pos)}:fill=#edc5d9,outline=#c63e83,radius=25"
    command.extend(["-highlightRegion", highlight_region_value])

    # --- Execute the command for this RNA ---
    print(command)
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
print(f"Output saved to '{OUTPUT_DIR}'.")
print("-------------------------------------------------")
