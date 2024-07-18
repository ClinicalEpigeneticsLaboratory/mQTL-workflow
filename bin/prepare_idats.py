#!/usr/bin/python3.10

import sys
import shutil
from os import makedirs
from os.path import join

from tqdm import tqdm
from glob import glob
import pandas as pd 

def validate(input_dir: str, sample_sheet: pd.DataFrame, array_position: str) -> None:
    expected_files_prefix = sample_sheet[array_position].unique()
    
    for prefix in tqdm(expected_files_prefix):
        grnIDAT = glob(join(input_dir, f"{prefix}*Grn.idat"))
        redIDAT = glob(join(input_dir, f"{prefix}*Red.idat"))

        if not grnIDAT:
            raise Exception(f"Grn IDAT for {prefix} defined in sample sheet column {array_position} not exists in {input_dir}")

        if len(grnIDAT) != 1:
            raise Exception(f"Duplicated Grn IDAT for {prefix} defined in sample sheet column {array_position}")
        
        if not redIDAT:
            raise Exception(f"Grn IDAT for {prefix} defined in sample sheet column {array_position} not exists in {input_dir}")

        if len(redIDAT) != 1:
            raise Exception(f"Duplicated Red IDAT for {prefix} defined in sample sheet column {array_position}")

    
def copy_idats(input_dir: str, sample_sheet: pd.DataFrame, array_position: str, dest: str) -> None:
    expected_files_prefix = sample_sheet[array_position].unique()
    
    for prefix in tqdm(expected_files_prefix):
        src = glob(join(input_dir, f"{prefix}*.idat"))
        
        if len(src) == 2:
            
            for file in src:
                shutil.copy(file, dest)

if __name__ == "__main__":
    input_dir, sample_sheet_path, array_position_col, dest = sys.argv[1:]
    
    print(f"""
    Processing IDATs dir: {input_dir}
    =================================
    Sample sheet: {sample_sheet_path}
    Array position column: {array_position_col}
    Destination directory: {dest}
    """)

    makedirs(dest)
    sample_sheet = pd.read_csv(sample_sheet_path)
    
    validate(input_dir, sample_sheet, array_position_col)
    copy_idats(input_dir, sample_sheet, array_position_col, dest)
