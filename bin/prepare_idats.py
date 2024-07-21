#!/usr/bin/python3.10

import sys
import shutil
from os import makedirs
from os.path import join

from tqdm import tqdm
from glob import glob
import pandas as pd


def copy_idats(
    input_dir: str, sample_sheet: pd.DataFrame, array_position: str, dest: str
) -> None:
    expected_files_prefix = sample_sheet[array_position].unique()

    for prefix in tqdm(expected_files_prefix):
        src = glob(join(input_dir, f"{prefix}*.idat"))

        if len(src) == 2:
            for file in src:
                shutil.copy(file, dest)


if __name__ == "__main__":
    input_dir, sample_sheet_path, array_position_col, dest = sys.argv[1:]

    print(
        f"""
    INPUT:
    =================================
    Processing IDATs dir: {input_dir}
    Sample sheet: {sample_sheet_path}

    PARAMS:
    =================================
    Array position column: {array_position_col}
    Destination directory: {dest}
    """
    )

    makedirs(dest)
    sample_sheet = pd.read_csv(sample_sheet_path)
    copy_idats(input_dir, sample_sheet, array_position_col, dest)
