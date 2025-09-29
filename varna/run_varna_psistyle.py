import pandas as pd
import os
import subprocess
import time

# replace all png with svg to produce vector based files

# --- Configuration ---
VARNA_JAR_PATH = os.path.expanduser('/home/users/rodell/varna/VARNAv3-93.jar')
PICKLE_FILE_PATH = 'consolidated_rna_data_withPUS7.pkl'

# Define the two output directories based on the PUS7 intersection
PUS7_OUTPUT_DIR = "varna_structures_PUS7intersection_png"
NON_PUS7_OUTPUT_DIR = "varna_structures_notPUS7intersection_png"
os.makedirs(PUS7_OUTPUT_DIR, exist_ok=True)
os.makedirs(NON_PUS7_OUTPUT_DIR, exist_ok=True)

print("--- VARNA High-Throughput Structure Generation ---")

# 1. Load the enriched DataFrame
try:
    df = pd.read_pickle(PICKLE_FILE_PATH)
    print(f"Successfully loaded enriched DataFrame from '{PICKLE_FILE_PATH}'.")
except FileNotFoundError:
    print(f"Error: The file '{PICKLE_FILE_PATH}' was not found. Please run the updated 'prepare_data.py' first.")
    exit()

total_rnas = len(df)
start_time = time.time()

# --- Main Loop: Iterate over every RNA in the DataFrame ---
for index, row in df.iterrows():
    rna_name = row['name']
    print(f"\nProcessing RNA {index + 1}/{total_rnas}: {rna_name}")

    # --- Set the correct output directory for this RNA ---
    has_pus7_site = pd.notna(row['psi_position'])
    if has_pus7_site:
        current_output_dir = PUS7_OUTPUT_DIR
    else:
        current_output_dir = NON_PUS7_OUTPUT_DIR

    # Prepare file paths
    safe_filename = rna_name.replace(',', '_').replace('/', '_')
    structure_filepath = os.path.join(current_output_dir, f"{safe_filename}.db")
    output_filepath = os.path.join(current_output_dir, f"{safe_filename}_structure.png")

    # Create the temporary .db file for this RNA
    with open(structure_filepath, "w") as f:
        f.write(f">{rna_name}\n")
        f.write(f"{row['varna_sequence']}\n") 
        f.write(f"{row['mea_structure']}\n")

    # --- Construct the VARNA command with your final theme ---
    command = [
        "java", "-cp", VARNA_JAR_PATH, "fr.orsay.lri.varna.applications.VARNAcmd",
        "-i", structure_filepath, "-o", output_filepath,
        
        # General SHAPE-based styles
        "-basesStyle1", "outline=#000000,label=#000000,fill=#e5e4e2",
        "-basesStyle2", "outline=#b26200,label=#b26200,fill=#ffc87c",
        "-basesStyle3", "outline=#6b1010,label=#6b1010,fill=#e0a3a3",
        
        # Combined styles for the Psi site
        "-basesStyle4", "fill=#dcdcdc,label=#c63e83,outline=#000000", # Low-reactivity Psi
        "-basesStyle5", "fill=#ffc87c,label=#c63e83,outline=#b26200", # Med-reactivity Psi
        "-basesStyle6", "fill=#e0a3a3,label=#c63e83,outline=#6b1010", # High-reactivity Psi

        # General options
        "-title", rna_name, "-resolution", "2.0", "-periodNum", "10"
    ]

    # --- Apply styles based on the data for this specific RNA ---
    low_bases = list(row['low_reactivity_bases'])
    med_bases = list(row['medium_reactivity_bases'])
    high_bases = list(row['high_reactivity_bases'])
    psi_pos = row['psi_position']

    # Handle the special Psi base first if it exists
    if has_pus7_site:
        print("  - Psi site found. Applying special layered highlighting.")
        psi_pos_str = str(int(psi_pos))

        if psi_pos_str in low_bases:
            command.extend(["-applyBasesStyle4on", psi_pos_str])
            low_bases.remove(psi_pos_str)
        elif psi_pos_str in med_bases:
            command.extend(["-applyBasesStyle5on", psi_pos_str])
            med_bases.remove(psi_pos_str)
        elif psi_pos_str in high_bases:
            command.extend(["-applyBasesStyle6on", psi_pos_str])
            high_bases.remove(psi_pos_str)
        
        # Draw the highlight box around the single Psi base
        region_str = f"{int(psi_pos)}-{int(psi_pos)}"
        command.extend(["-highlightRegion", f"{region_str}:fill=#edc5d9,outline=#c63e83,radius=25"])

    # Apply the general SHAPE styles to all remaining bases
    if low_bases: command.extend(["-applyBasesStyle1on", ",".join(low_bases)])
    if med_bases: command.extend(["-applyBasesStyle2on", ",".join(med_bases)])
    if high_bases: command.extend(["-applyBasesStyle3on", ",".join(high_bases)])

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
print(f"Output saved to '{PUS7_OUTPUT_DIR}' and '{NON_PUS7_OUTPUT_DIR}'.")
print("-------------------------------------------------")