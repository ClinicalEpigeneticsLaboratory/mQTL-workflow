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
    manifest = manifest[["Name", "Chr", "MapInfo"]]
    
    manifest = manifest.dropna().set_index("Name")
    manifest = manifest.loc[
        list(set(manifest.index).intersection(set(genomic_frame.index)))
    ]

    print(manifest)
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

    manifest = pd.read_csv(manifest, skiprows=7, low_memory=False).set_index("Name")

    manifest = manifest.loc[
        list(set(manifest.index).intersection(set(methylation_frame.index)))
    ]

    manifest = manifest[["CHR", "MAPINFO"]]
    manifest["MAPINFO"] = manifest["MAPINFO"].astype(int)

    for chr_ in manifest.CHR.unique():
        partial = manifest[manifest.CHR == chr_]
        map[f"chr{chr_}"] = partial

    return map


def encode(values: pd.Series) -> dict[str, int]:
    order = sorted(values.unique())
    numerical = list(range(len(order)))

    return dict(zip(order, numerical))


def analyse(chr_):
    stats = []
    print(methylation_map.keys())
    
    snps_positions_per_chr = gsa_map[chr_]
    cpgs_positions_per_chr = methylation_map[chr_]

    for snpID in tqdm(snps_positions_per_chr.index):
        snp_pos = snps_positions_per_chr.loc[snpID, "MAPINFO"]
        snp_data = genotype_frame.loc[snpID]

        mapper = encode(snp_data)
        snp_data = snp_data.map(mapper)

        cpgs_in_range = cpgs_positions_per_chr["MAPINFO"] - snp_pos
        cpgs_in_range = cpgs_in_range[
            cpgs_in_range.abs() <= int(distance)
        ].index.tolist()

        if not cpgs_in_range:
            continue

        cpgs_in_range = methylation_frame.loc[cpgs_in_range]

        for cpg, cpg_data in cpgs_in_range.iterrows():
            cpg_pos = cpgs_positions_per_chr.loc[cpg, "MAPINFO"]

            model = sts.linregress(x=snp_data, y=cpg_data)
            slope, r, pval = model.slope, model.rvalue, model.pvalue
            stats.append(
                {
                    "CHR": chr_,
                    "snp": snpID,
                    "snp POS": snp_pos,
                    "cpg": cpg,
                    "cpg POS": cpg_pos,
                    "slope": abs(slope),
                    "R2": r**2,
                    "p-value": pval
                }
            )

    data = pd.DataFrame(stats)
    _, data["FDR"] = fdrcorrection(data["p-value"], method="n")
    
    data["distance"] = data["snp POS"] - data["cpg POS"]
    data.distance = data.distance.abs()

    return data


if __name__ == "__main__":
    (
        genotype_frame,
        methylation_frame,
        gsa_manifest,
        methylation_manifest,
        conversion_file,
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
    Conversion file: {conversion_file}

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
    
    conversion_file = pd.read_table(conversion_file).set_index("Name")
    results = pd.merge(conversion_file, results, how="right", left_index=True, right_on="snp")
    
    results.to_parquet("mQTL.parquet")
