#!/usr/bin/python3
import sys

import pandas as pd


def recode_sample_names(
    df: pd.DataFrame,
    sample_sheet: pd.DataFrame,
    array_position_col: str,
    sample_name_col: str,
) -> pd.DataFrame:
    common = list(set(sample_sheet[array_position_col]).intersection(set(df.columns)))
    df = df[common]

    mapper = dict(zip(sample_sheet[array_position_col], sample_sheet[sample_name_col]))
    df.columns = df.columns.map(lambda x: mapper[x] if x in mapper.keys() else x)
    return df


if __name__ == "__main__":
    mynorm_path, sample_sheet_path, array_position_col, sample_name_col = sys.argv[1:]

    print(
        f"""
    Processing mynorm file: {mynorm_path}
    =================================
    Sample sheet: {sample_sheet_path}
    Array position column: {array_position_col}
    Sample name column: {sample_name_col}
    """
    )

    mynorm = pd.read_parquet(mynorm_path)

    if "CpG" in mynorm.columns:
        mynorm = mynorm.set_index("CpG")

    sample_sheet = pd.read_csv(sample_sheet_path)

    mynorm = recode_sample_names(
        mynorm, sample_sheet, array_position_col, sample_name_col
    )
    mynorm.dropna().to_parquet(f"mynorm.parquet")
