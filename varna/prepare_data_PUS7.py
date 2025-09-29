import pandas as pd
import os
import re

# --- Configuration ---
FOLD_DIR = '/scratch/groups/nicolemm/rodell/SHAPE/InVitro_2repHighQ_Analysis/RNAfold_output'
SHAPE_DIR = '/scratch/groups/nicolemm/rodell/SHAPE/InVitro_2repHighQ_Analysis/averaged_shape_files'
PUS7_CSV_PATH = '/scratch/groups/nicolemm/rodell/SHAPE/InVitro_2repHighQ_Analysis/varna/PUS7sites_DiffLengths.csv'

# --- Data Loading and Parsing ---
all_rna_data = []
print(f"Starting to process files in {FOLD_DIR}...")

# (This part is the same as before, loading the base .fold and .shape data)
for fold_filename in os.listdir(FOLD_DIR):
    if not fold_filename.endswith('.fold'):
        continue
    fold_filepath = os.path.join(FOLD_DIR, fold_filename)
    try:
        with open(fold_filepath, 'r') as f:
            lines = [line.strip() for line in f.readlines() if line.strip()]
        header = lines[0]
        if not header.startswith('>'): continue
        full_rna_name = header[1:]
        sequence = lines[1]
        mea_line = lines[5]
        if 'MEA=' not in mea_line: continue
        mea_structure = mea_line.split(' ')[0]
        shape_file_base_name = full_rna_name.replace(',', '_')
        shape_filename = f"{shape_file_base_name}.shape"
        shape_filepath = os.path.join(SHAPE_DIR, shape_filename)
        if not os.path.exists(shape_filepath): continue
        shape_df = pd.read_csv(shape_filepath, sep='\s+', header=None, names=['position', 'reactivity'])
        shape_reactivities = shape_df['reactivity'].tolist()
        if len(sequence) != len(shape_reactivities): continue
        
        all_rna_data.append({
            'name': full_rna_name,
            'sequence': sequence,
            'mea_structure': mea_structure,
            'shape_reactivities': shape_reactivities
        })
    except Exception as e:
        print(f"  [Error] Failed to process file {fold_filename}: {e}")

df = pd.DataFrame(all_rna_data)
print(f"\nLoaded {len(df)} initial RNA structures.")

# --- Load and merge PUS7 pseudouridine site data ---
print(f"Loading PUS7 site data from {PUS7_CSV_PATH}...")
pus7_df = pd.read_csv(PUS7_CSV_PATH)

# Sanitize the 'name' column in BOTH DataFrames to create a common merge key.
df['merge_key'] = df['name'].str.replace(',', '_')
pus7_df['merge_key'] = pus7_df['name'].str.replace(',', '_')

# Keep only the columns we need from the PUS7 data for the merge
pus7_to_merge = pus7_df[['merge_key', 'original_psipos']]

# Perform a left merge using the new, common 'merge_key'.
df = pd.merge(df, pus7_to_merge, on='merge_key', how='left')

# Clean up the temporary merge key column
df = df.drop(columns=['merge_key'])
print(f"Merged PUS7 data. Found {df['original_psipos'].notna().sum()} RNAs with a known Psi site.")


# --- Pre-process data for VARNA visualization ---
print("Pre-processing data for VARNA...")

# Initialize new columns
df['varna_sequence'] = df['sequence']
df['psi_position'] = None
df['highlight_region'] = None
df['low_reactivity_bases'] = [[] for _ in range(len(df))]
df['medium_reactivity_bases'] = [[] for _ in range(len(df))]
df['high_reactivity_bases'] = [[] for _ in range(len(df))]

for index, row in df.iterrows():
    # 1. Process Psi site information, if it exists
    if pd.notna(row['original_psipos']):
        psipos = int(row['original_psipos'])
        df.at[index, 'psi_position'] = psipos
        
        # Define the 5-mer region for the highlight box (Psi +/- 2)
        start_pos = max(1, psipos - 2)
        end_pos = min(len(row['sequence']), psipos + 2)
        df.at[index, 'highlight_region'] = f"{start_pos}-{end_pos}"
        
        # Create the VARNA sequence with 'P' for Pseudouridine
        seq_list = list(row['sequence'])
        if seq_list[psipos - 1] == 'U': # Sanity check
             seq_list[psipos - 1] = 'P'
        df.at[index, 'varna_sequence'] = "".join(seq_list)

    # 2. Process SHAPE reactivity for all RNAs
    low_bases, medium_bases, high_bases = [], [], []
    for i, reactivity in enumerate(row['shape_reactivities']):
        base_number = i + 1
        if reactivity < 0.3:
            low_bases.append(str(base_number))
        elif reactivity <= 0.7:
            medium_bases.append(str(base_number))
        else:
            high_bases.append(str(base_number))
            
    df.at[index, 'low_reactivity_bases'] = low_bases
    df.at[index, 'medium_reactivity_bases'] = medium_bases
    df.at[index, 'high_reactivity_bases'] = high_bases

# --- Save the enriched DataFrame ---
output_pickle_file = "consolidated_rna_data_withPUS7.pkl"
df.to_pickle(output_pickle_file)
print(f"\nEnriched DataFrame with PUS7 and SHAPE categories saved to '{output_pickle_file}'")