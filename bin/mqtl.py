#!/usr/bin/python3.10

import sys
from functools import partial
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


def encode(values: pd.Series) -> dict[str, int]:
    order = sorted(values.unique())
    numerical = list(range(len(order)))

    return dict(zip(order, numerical))


def analyse(chr_: str, 
            min_distance: int,
            max_distance: int, 
            gsa_map: dict[str, pd.DataFrame], 
            methylation_map: dict[str, pd.DataFrame], 
            genotype_frame: pd.DataFrame, 
            methylation_frame: pd.DataFrame, 
            sample_sheet: pd.DataFrame, 
            covariates: list[str]
           ) -> pd.DataFrame:

    stats = []
    snps_positions_per_chr = gsa_map[chr_]
    cpgs_positions_per_chr = methylation_map[chr_]

    for snpID in tqdm(snps_positions_per_chr.index[:100]):
        snp_pos = snps_positions_per_chr.loc[snpID, "MAPINFO"]
        snp_data = genotype_frame.loc[snpID]
        
        mapper = encode(snp_data)
        snp_data = snp_data.map(mapper)

        cpgs_in_range = cpgs_positions_per_chr.loc[
            (cpgs_positions_per_chr["MAPINFO"] - snp_pos).abs().between(int(min_distance), int(max_distance))
        ]

        if cpgs_in_range.empty:
            continue

        cpg_data_in_range = methylation_frame.loc[cpgs_in_range.index]
        exog_data = pd.concat((sample_sheet, snp_data), axis=1, join="inner")

        for cpg in cpg_data_in_range.index:
            cpg_data = cpg_data_in_range.loc[cpg]
            model = sms.RLM(cpg_data, exog_data, M=sms.robust.norms.HuberT()).fit()
            
            slope = model.params.get(snpID)
            pval = model.pvalues.get(snpID)

            stats.append(
                {
                    "CHR": chr_,
                    "snp": snpID,
                    "snp POS": snp_pos,
                    "cpg": cpg,
                    "cpg POS": cpgs_positions_per_chr.loc[cpg, "MAPINFO"],
                    "|slope|": abs(slope),
                    "p-value": pval,
                    "n": len(common_samples),
                    "Model": f"{cpg} ~ {' + '.join(exog_data.columns)}"
                }
            )

    data = pd.DataFrame(stats)
    _, data["FDR"] = fdrcorrection(data["p-value"], method="n")
    data["distance"] = (data["snp POS"] - data["cpg POS"]).abs()

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

    methylation_map = generate_methylation_map(methylation_manifest, methylation_frame)
    gsa_map = generate_gsa_map(gsa_manifest, genotype_frame)

    common_samples = list(set.intersection(set(methylation_frame.columns), set(genotype_frame.columns)))
    genotype_frame = genotype_frame[common_samples]
    methylation_frame = methylation_frame[common_samples]

    sample_sheet = pd.read_csv(sample_sheet).set_index(sample_name)
    sample_sheet = sample_sheet.loc[common_samples]
    sample_sheet.insert(loc=0, column="Intercept", value=1)

    if covariates:
        covariates = [covar.strip() for covar in covariates.split(",")]
        sample_sheet = sample_sheet[["Intercept", *covariates]]

    else:
        sample_sheet = sample_sheet[["Intercept"]]
    
    f = partial(analyse, 
                      min_distance=min_distance, 
                      max_distance=max_distance, 
                      gsa_map=gsa_map, 
                      methylation_map=methylation_map, 
                      genotype_frame=genotype_frame, 
                      methylation_frame=methylation_frame, 
                      sample_sheet=sample_sheet, 
                      covariates=covariates
                     )
    
    with Pool(int(nCPU)) as p:
        results = p.map(f, gsa_map.keys())

    results = pd.concat(results)
    
    conversion_file = pd.read_table(conversion_file).set_index("Name")
    results = pd.merge(conversion_file, results, how="right", left_index=True, right_on="snp")
    
    results.to_parquet("mQTL.parquet")
