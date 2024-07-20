#!/usr/bin/python3.10

import sys

import pandas as pd
from tqdm import tqdm
import numpy as np


def preprocess_traw_file(traw_path: str) -> pd.DataFrame:
    df = pd.read_table(traw_path)
    df = df.drop("(C)M", axis=1)
    df = df.dropna()

    return df


def preprocess_sample_sheet(
    sample_sheet_path: str, array_position_col: str, sample_name_col: str
) -> pd.DataFrame:
    sample_sheet = pd.read_csv(sample_sheet_path)

    if not array_position_col in sample_sheet.columns:
        raise Exception(f"{array_position_col} is not valid sample sheet column")

    if not sample_name_col in sample_sheet.columns:
        raise Exception(f"{sample_name_col} is not valid sample sheet column")

    sample_sheet = sample_sheet[[array_position_col, sample_name_col]].dropna()

    return sample_sheet


def recode_sample_names(
    df: pd.DataFrame,
    sample_sheet: pd.DataFrame,
    array_position_col: str,
    sample_name_col: str,
) -> pd.DataFrame:
    common = list(set(sample_sheet[array_position_col]).intersection(set(df.columns)))

    df = df[["SNP", "COUNTED", "ALT", *common]]

    mapper = dict(zip(sample_sheet[array_position_col], sample_sheet[sample_name_col]))
    df.columns = df.columns.map(lambda x: mapper[x] if x in mapper.keys() else x)
    return df


def recode_notation(df: pd.DataFrame):
    encoded = []
    for id_, row in tqdm(df.iterrows(), total=df.shape[0]):
        a, b = row["COUNTED"], row["ALT"]
        mapper = {0.0: f"{b}/{b}", 1.0: f"{a}/{b}", 2.0: f"{a}/{a}"}

        encoded.append(row.replace(mapper).to_frame().T)

    return pd.concat(encoded)


if __name__ == "__main__":
    traw_path, sample_sheet_path, array_position_col, sample_name_col = sys.argv[1:]

    print(
        f"""
    INPUT:
    =================================
    Processing traw file: {traw_path}
    Sample sheet: {sample_sheet_path}

    PARAMS:
    =================================
    Array position column: {array_position_col}
    Sample name column: {sample_name_col}
    """
    )

    traw = preprocess_traw_file(traw_path)
    sample_sheet = preprocess_sample_sheet(
        sample_sheet_path, array_position_col, sample_name_col
    )

    traw = recode_sample_names(traw, sample_sheet, array_position_col, sample_name_col)
    traw = recode_notation(traw)
    traw = traw.set_index("SNP")

    traw = traw.drop(["COUNTED", "ALT"], axis=1)
    traw = traw.loc[traw.nunique(axis=1) > 1,]

    traw.to_parquet("genotype_table.parquet")
