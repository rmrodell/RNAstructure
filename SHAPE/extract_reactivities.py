#!/usr/bin/env python3
import pandas as pd
import argparse
from pathlib import Path
import sys

def concatenate_shape_files(input_dir, output_file):
    """
    Combines SHAPE reactivity data from .shape files in a directory
    into a single CSV file.

    Args:
        input_dir (str): The path to the directory containing .shape files.
        output_file (str): The path for the output CSV file.
    """
    input_path = Path(input_dir)
    if not input_path.is_dir():
        print(f"Error: Input directory not found at '{input_dir}'", file=sys.stderr)
        sys.exit(1)

    shape_files = list(input_path.glob('*.shape'))

    if not shape_files:
        print(f"Warning: No .shape files found in '{input_dir}'", file=sys.stderr)
        return

    print(f"Found {len(shape_files)} .shape files to process...")

    all_dataframes = []

    for file in shape_files:
        file_basename = file.stem
        print(f"  - Processing {file.name}...")

        try:
            df = pd.read_csv(
                file,
                sep='\s+',
                header=None,
                names=['position', 'reactivity'],
                comment='#' # Ignore any lines that start with #
            )
            
            df['name'] = file_basename
            
            all_dataframes.append(df)

        except pd.errors.EmptyDataError:
            print(f"    Warning: {file.name} is empty and will be skipped.")
        except Exception as e:
            print(f"    Error processing {file.name}: {e}")

    if not all_dataframes:
        print("No valid data was processed. Output file will not be created.")
        return

    # Concatenate all DataFrames in the list into a single one
    long_df = pd.concat(all_dataframes, ignore_index=True)
    
    try:
        pivoted_df = long_df.pivot(index='name', columns='position', values='reactivity')
    except Exception as e:
        print(f"Error: Could not pivot data. This can happen if there are duplicate "
              f"position entries within a single file. Details: {e}", file=sys.stderr)
        sys.exit(1)

    # Remove the name of the columns index ('position')
    pivoted_df.columns.name = None
    # Move the index ('name') to be a regular column
    pivoted_df = pivoted_df.reset_index()
    
    pivoted_df.to_csv(output_file, index=False)
    
    print(f"\nSuccess! SHAPE reactivities saved to '{output_file}'")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Concatenate SHAPE reactivity data from multiple .shape files into a single CSV.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "input_dir",
        help="Directory containing the .shape files."
    )
    parser.add_argument(
        "output_file",
        help="Path for the output .csv file."
    )
    
    args = parser.parse_args()
    
    concatenate_shape_files(args.input_dir, args.output_file)