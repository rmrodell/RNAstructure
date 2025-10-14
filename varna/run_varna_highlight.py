import pandas as pd
import os
import subprocess
import time

# replace all png with svg to produce vector based files

# --- Configuration ---
VARNA_JAR_PATH = os.path.expanduser('/home/users/rodell/varna/VARNAv3-93.jar')

INPUT_FILE_PATH = '/scratch/users/rodell/20251014_pool2/RNAfold_fold_summary.csv'  # <--- UPDATE THIS TO YOUR FILENAME
HIGHLIGHT_POSITION = 66

OUTPUT_DIR = "varna_structures_svg"
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("--- VARNA High-Throughput Structure Generation ---")

# 1. Load the DataFrame
try:
    df = pd.read_csv(INPUT_FILE_PATH) 
    print(f"Successfully loaded CSV file from '{INPUT_FILE_PATH}'.")
except FileNotFoundError:
    print(f"Error: The file '{INPUT_FILE_PATH}' was not found.")
    exit()

total_rnas = len(df)
start_time = time.time()

# --- Main Loop: Iterate over every RNA in the DataFrame ---
for index, row in df.iterrows():
    rna_name = row['sequence_id']
    print(f"\nProcessing RNA {index + 1}/{total_rnas}: {rna_name}")

    # Prepare file paths
    safe_filename = rna_name.replace(',', '_').replace('/', '_')
    structure_filepath = os.path.join(OUTPUT_DIR, f"{safe_filename}.db")
    output_filepath = os.path.join(OUTPUT_DIR, f"{safe_filename}.svg")

    # Create the temporary .db file for this RNA
    with open(structure_filepath, "w") as f:
        f.write(f">{rna_name}\n")
        f.write(f"{row['sequence']}\n") 
        f.write(f"{row['mea_structure']}\n")

    highlight_region_value = f"{HIGHLIGHT_POSITION}-{HIGHLIGHT_POSITION}:fill=#edc5d9,outline=#c63e83,radius=25"

    # --- Construct the VARNA command with your final theme ---
    command = [
        "java", "-cp", VARNA_JAR_PATH, "fr.orsay.lri.varna.applications.VARNAcmd",
        "-i", structure_filepath, 
        "-o", output_filepath,
        
        # highlight residue style
        "-basesStyle2", "fill=#dcdcdc,label=#c63e83,outline=#000000",
        "-applyBasesStyle2on", str(HIGHLIGHT_POSITION),

        # default style
        "-basesStyle1", "outline=#000000,label=#000000,fill=#e5e4e2",
        "-applyBasesStyle1on", "1-1000",

        # General options
        "-title", rna_name, "-resolution", "2.0", "-periodNum", "10",

        "-highlightRegion", highlight_region_value
    ]

    # --- Execute the command for this RNA ---
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