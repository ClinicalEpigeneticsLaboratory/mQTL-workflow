#!/usr/bin/python3.10

import sys
from multiprocessing import Pool

import pandas as pd
from tqdm import tqdm
import statsmodels.api as sms
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


def generate_interaction_map(
    gsa_map: dict, methylation_map: dict, min_distance: int, max_distance: int
) -> dict:
    map = {}
    
    for chr in gsa_map.keys():
        map_per_chr = []
        snps_positions_per_chr = gsa_map[chr]
        cpgs_positions_per_chr = methylation_map[chr]

        for snpID in tqdm(snps_positions_per_chr.index, desc=chr):
            snp_pos = snps_positions_per_chr.loc[snpID, "MAPINFO"]
            cpgs_in_range = cpgs_positions_per_chr.loc[
                (cpgs_positions_per_chr["MAPINFO"] - snp_pos)
                .abs()
                .between(int(min_distance), int(max_distance))
            ]

            if cpgs_in_range.empty:
                continue

            for cpgID in cpgs_in_range.index:
                map_per_chr.append([snpID, cpgID])

        map[chr] = map_per_chr

    return map


def encode(values: pd.Series) -> dict[str, int]:
    order = sorted(values.unique())
    numerical = list(range(len(order)))

    mapper = dict(zip(order, numerical))
    return values.map(mapper)


def annotate(
    mqtls: pd.DataFrame,
    conversion_file: str,
    gsa_manifest: str,
    methylation_manifest: str,
) -> pd.DataFrame:
    conversion_file = pd.read_table(conversion_file).set_index("Name")
    gsa_manifest = pd.read_csv(gsa_manifest, skiprows=7, low_memory=False).set_index(
        "Name"
    )
    methylation_manifest = pd.read_csv(
        methylation_manifest, skiprows=7, low_memory=False
    ).set_index("Name")

    mqtls = pd.merge(
        mqtls, conversion_file, how="left", left_on="snp", right_index=True
    )
    mqtls = pd.merge(
        mqtls,
        gsa_manifest.rename({"MapInfo": "snp pos"}, axis=1)["snp pos"],
        how="left",
        left_on="snp",
        right_index=True,
    )
    mqtls = pd.merge(
        mqtls,
        methylation_manifest.rename({"MAPINFO": "cpg pos"}, axis=1)["cpg pos"],
        how="left",
        left_on="cpg",
        right_index=True,
    )

    mqtls["|distance|"] = (mqtls["snp pos"] - mqtls["cpg pos"]).abs()
    return mqtls


def prepare_chunks(interactions_map: dict, genotype_frame: pd.DataFrame, methylation_frame: pd.DataFrame, sample_sheet) -> list:
    chunks = []
    for chr_ in [f"chr{i+1}" for i in range(22)]:
        interactions = interactions_map[chr_]

        rss = list({x[0] for x in interactions})
        cpgs = list({x[1] for x in interactions})
        
        chunks.append([chr_, interactions, genotype_frame.loc[rss], methylation_frame.loc[cpgs], sample_sheet])

    return chunks


def analyse(data: list) -> pd.DataFrame:
    chr_, interactions, genotype_frame, methylation_frame, sample_sheet = data
    stats = []
    
    for snpID, cpgID in tqdm(interactions):
        exog = pd.concat((genotype_frame.loc[snpID], sample_sheet), axis=1)
        model = sms.RLM(methylation_frame.loc[cpgID], exog, M=sms.robust.norms.HuberT()).fit()

        slope = model.params[snpID]
        pval = model.pvalues[snpID]

        stats.append(
            {
                "CHR": chr_,
                "snp": snpID,
                "cpg": cpgID,
                "|slope|": abs(slope),
                "p-value": pval,
                "Model": f"{cpgID} ~ {' + '.join([*sample_sheet.columns.tolist(), snpID])}",
            }
        )

    data = pd.DataFrame(stats)
    _, data["FDR"] = fdrcorrection(data["p-value"], method="n")
    return data


if __name__ == "__main__":
    (
        genotype_frame,
        methylation_frame,
        gsa_manifest,
        methylation_manifest,
        conversion_file,
        sample_sheet,
        sample_name,
        covariates,
        min_distance,
        max_distance,
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
    Sample sheet file: {sample_sheet}
    Sample name column: {sample_name}
    Co-variates: {covariates}
    
    OPTIONS:
    =================================
    Min distance threshold [bp]: {min_distance}
    Max distance threshold [bp]: {max_distance}
    CPUs: {nCPU}
    """
    )
    methylation_frame = pd.read_parquet(methylation_frame)
    if "CpG" in methylation_frame.columns:
        methylation_frame = methylation_frame.set_index("CpG")

    genotype_frame = pd.read_parquet(genotype_frame)
    if "SNP" in genotype_frame.columns:
        genotype_frame = genotype_frame.set_index("SNP")
    genotype_frame = genotype_frame.apply(encode, axis=1)

    methylation_map = generate_methylation_map(methylation_manifest, methylation_frame)
    gsa_map = generate_gsa_map(gsa_manifest, genotype_frame)
    interactions_map = generate_interaction_map(
        gsa_map, methylation_map, min_distance, max_distance
    )

    sample_sheet = pd.read_csv(sample_sheet).set_index(sample_name)
    common_samples = list(
        set.intersection(set(methylation_frame.columns), set(genotype_frame.columns), set(sample_sheet.index))
    )
    
    genotype_frame = genotype_frame[common_samples]
    methylation_frame = methylation_frame[common_samples]
    sample_sheet = sample_sheet.loc[common_samples]
    sample_sheet.insert(loc=0, column="Intercept", value=1)

    if covariates:
        covariates = [covar.strip() for covar in covariates.split(",")]
        sample_sheet = sample_sheet[["Intercept", *covariates]]

    else:
        sample_sheet = sample_sheet[["Intercept"]]

    chunks = prepare_chunks(interactions_map, genotype_frame, methylation_frame, sample_sheet)
    with Pool(int(nCPU)) as p:
        results = p.map(analyse, chunks)

    results = pd.concat(results)
    results = results.reset_index(drop=True)
    
    results = annotate(results, conversion_file, gsa_manifest, methylation_manifest)
    results.to_parquet("mQTL.parquet")
