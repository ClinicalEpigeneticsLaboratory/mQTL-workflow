#!/usr/bin/python3.10

import sys
from functools import partial
from multiprocessing import Pool

import pandas as pd
from tqdm import tqdm
import scipy.stats as sts
from statsmodels.stats.multitest import fdrcorrection


def generate_gsa_map(
    manifest: str, genomic_frame: pd.DataFrame
) -> dict[str, pd.DataFrame]:
    map = {}

    manifest = pd.read_csv(manifest, skiprows=7, low_memory=False)

    manifest["rsID"] = manifest.IlmnID.str.extract(r"(rs\d+)")
    manifest = manifest[["rsID", "Chr", "MapInfo"]]

    manifest = manifest.dropna().set_index("rsID")
    manifest = manifest.loc[
        list(set(manifest.index).intersection(set(genomic_frame.index)))
    ]

    manifest.columns = ["CHR", "MAPINFO"]
    manifest["MAPINFO"] = manifest["MAPINFO"].astype(int)

    for chr_ in manifest.CHR.unique():
        partial = manifest[manifest.CHR == chr_]
        map[f"chr{chr_}"] = partial

    return map


def generate_methylation_map(
    manifest: str, methylation_frame: pd.DataFrame
) -> dict[str, pd.DataFrame]:
    map = {}

    manifest = pd.read_csv(manifest, skiprows=7, low_memory=False, index_col=0)
    manifest = manifest.loc[
        list(set(manifest.index).intersection(set(methylation_frame.index)))
    ]

    manifest = manifest[["CHR", "MAPINFO"]]
    manifest["MAPINFO"] = manifest["MAPINFO"].astype(int)

    for chr_ in manifest.CHR.unique():
        partial = manifest[manifest.CHR == chr_]
        map[f"chr{chr_}"] = partial

    return map


def encode(values: pd.Series):
    order = sorted(values.unique())
    numerical = list(range(len(order)))

    return dict(zip(order, numerical))


def analyse(chr_):
    stats = []

    rss_positions_per_chr = gsa_map[chr_]
    cpgs_positions_per_chr = methylation_map[chr_]

    for rsID in tqdm(rss_positions_per_chr.index):
        rs_pos = rss_positions_per_chr.loc[rsID, "MAPINFO"]
        rs_data = genotype_frame.loc[rsID]

        mapper = encode(rs_data)
        rs_data = rs_data.replace(mapper)

        cpgs_in_range = cpgs_positions_per_chr["MAPINFO"] - rs_pos
        cpgs_in_range = cpgs_in_range[
            cpgs_in_range.abs() <= int(distance)
        ].index.tolist()

        if not cpgs_in_range:
            continue

        cpgs_in_range = methylation_frame.loc[cpgs_in_range]

        for cpg, cpg_data in cpgs_in_range.iterrows():
            cpg_pos = cpgs_positions_per_chr.loc[cpg, "MAPINFO"]

            model = sts.linregress(x=rs_data, y=cpg_data)
            slope, r, pval = model.slope, model.rvalue, model.pvalue
            stats.append(
                {
                    "CHR": chr_,
                    "rs": rsID,
                    "rs POS": rs_pos,
                    "cpg": cpg,
                    "cpg POS": cpg_pos,
                    "slope": abs(slope),
                    "R2": r**2,
                    "p-value": pval,
                }
            )

    data = pd.DataFrame(stats)
    _, data["FDR"] = fdrcorrection(data["p-value"], method="n")
    data["distance"] = data["rs POS"] - data["cpg POS"]
    data.distance = data.distance.abs()

    return data


if __name__ == "__main__":
    (
        genotype_frame,
        methylation_frame,
        gsa_manifest,
        methylation_manifest,
        distance,
        nCPU,
    ) = sys.argv[1:]

    print(
        f"""
    INPUT:
    =================================
    Processing genotype frame file: {genotype_frame}
    Processing methylation frame file: {methylation_frame}
    GSA manifest: {gsa_manifest}
    Methylation manifest: {methylation_manifest}

    OPTIONS:
    =================================
    Distance threshold [bp]: {distance}
    CPUs: {nCPU}
    """
    )

    methylation_frame = pd.read_parquet(methylation_frame)
    if "CpG" in methylation_frame.columns:
        methylation_frame = methylation_frame.set_index("CpG")

    genotype_frame = pd.read_parquet(genotype_frame)
    if "SNP" in genotype_frame.columns:
        genotype_frame = genotype_frame.set_index("SNP")

    methylation_map = generate_methylation_map(methylation_manifest, methylation_frame)
    gsa_map = generate_gsa_map(gsa_manifest, genotype_frame)

    common_samples = list(set.intersection(set(methylation_frame.columns), set(genotype_frame.columns)))
    genotype_frame = genotype_frame[common_samples]
    methylation_frame = methylation_frame[common_samples]
    
    with Pool(int(nCPU)) as p:
        results = p.map(analyse, gsa_map.keys())

    results = pd.concat(results)
    results.to_parquet("mQTL.parquet")
