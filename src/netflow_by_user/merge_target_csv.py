#!/usr/bin/env python3
import os
import sys
import argparse
import glob
import numpy as np
from astropy.table import Table, vstack

def main():
    parser = argparse.ArgumentParser(description="Merge all ECSV files in a directory, drop duplicates, and save to the parent directory.")
    parser.add_argument("directory", help="Path to the directory containing ECSV files to merge")
    args = parser.parse_args()

    target_dir = os.path.normpath(args.directory)
    if not os.path.isdir(target_dir):
        print(f"Error: Directory '{target_dir}' does not exist.")
        sys.exit(1)

    # Find all ECSV files
    ecsv_files = glob.glob(os.path.join(target_dir, "*.ecsv"))
    if not ecsv_files:
        print(f"No ECSV files found in '{target_dir}'.")
        sys.exit(0)

    print(f"Found {len(ecsv_files)} ECSV files in '{target_dir}'. Reading...")

    tables = []
    for f in ecsv_files:
        try:
            tbl = Table.read(f, format="ascii.ecsv")
            tables.append(tbl)
        except Exception as e:
            print(f"  Warning: Failed to read {f}: {e}")

    if not tables:
        print("No tables were successfully loaded.")
        sys.exit(1)

    print("Merging tables...")
    try:
        merged_table = vstack(tables)
    except Exception as e:
        print(f"Error: Failed to vstack tables: {e}")
        print("This usually happens if columns or datatypes do not match across files.")
        sys.exit(1)

    # Identify primary key for deduplication
    id_col = None
    for col in ['obj_id', 'source_id', 'fluxstd_id', 'sky_id', 'ob_code']:
        if col in merged_table.colnames:
            id_col = col
            break

    print(f"Deduplicating targets (using key column: '{id_col if id_col else 'None'}')...")
    initial_rows = len(merged_table)
    
    if id_col is not None:
        # Deduplicate based on ID column, preserving first occurrence
        _, unique_indices = np.unique(merged_table[id_col].data, return_index=True)
        unique_indices.sort()
        deduped_table = merged_table[unique_indices]
    else:
        # Fallback to pandas for row-wise exact match deduplication
        df = merged_table.to_pandas()
        df = df.drop_duplicates()
        deduped_table = Table.from_pandas(df)

    final_rows = len(deduped_table)
    print(f"Merged {initial_rows} rows down to {final_rows} unique rows (removed {initial_rows - final_rows} duplicates).")

    # Determine output path
    parent_dir = os.path.dirname(os.path.abspath(target_dir))
    dir_name = os.path.basename(target_dir)
    output_path = os.path.join(parent_dir, f"{dir_name}.ecsv")

    print(f"Writing merged ECSV to: {output_path}...")
    deduped_table.write(output_path, format="ascii.ecsv", overwrite=True)
    print("Done!")

if __name__ == "__main__":
    main()
