# mQTL-workflow
This Nextflow project is designed to identify methylation quantitative trait loci (mQTL), which are interactions between the genome and the methylome. The project consists of three main workflows:

1. `flow_one` - Preprocess SNP arrays data (Illumina).
2. `flow_two` - Preprocess methylation arrays data (Illumina).
3. `flow_three` - Identify mQTLs.

### Getting Started

#### Prerequisites
Ensure you have the following software installed:

- Java >11, <=22
- Nextflow
- Docker

#### Docker Container
This projects contains docker image mainly for Python3.10 as well as R4.4 used as execution environment for Nextflow pipelines.
To build image, use the following command:

```
docker build -t mqtl-image .
```

To validate builds:

```
docker run -it mqtl-image

EXPECTED OUTPUT:

Running container using following dependencies:
Python: Python 3.10.0
R: R version 4.4.1 (2024-06-14) -- "Race for Your Life"
bgzip: bgzip (htslib) 1.9 Copyright (C) 2018 Genome Research Ltd.
tabix: tabix (htslib) 1.9 Copyright (C) 2018 Genome Research Ltd.
BCFtools: bcftools 1.9 Using htslib 1.9 Copyright (C) 2018 Genome Research Ltd. License Expat: The MIT/Expat license This is free software: you are free to change and redistribute it. There is NO WARRANTY, to the extent permitted by law.
Illumina CLI: Array Analysis CLI 2.1.0

[...]
```

### Workflow one: SNP Array Data Preprocessing
The first workflow preprocesses SNP array data from Illumina microarrays. It involves the following steps:

```
0. Sanity check
1. Convert idat to gtc: Convert raw idat files to gtc format.
2. Convert gtc to vcf: Convert gtc files to VCF format.
3. Index and merge VCF files: Index and merge VCF files into a BCF file.
4. Filtering and formatting:
    1. Use PLINK to filter the merged BCF file.
    2. Apply filters: minor allele frequency (--maf 0.05), Hardy-Weinberg equilibrium (--hwe 1e-50), genotype missingness (--geno 0.1), and autosome restriction.
    3. Extract only bi-allelic loci.
    4. Recode to a transposed format.
5. Exporting results files (*.bcf, *.frq, *.traw as well as recoded traw file - *.parquet).
```

Workflow basic command:

```
nextflow run flow_one.nf --reference_fa <path> --bmp_manifest <path> --csv_manifest <path> --cluster_file <path> --gsa_idats_dir <path> --sample_sheet <path> --array_position <str> --sample_name <str> --results_dir <path>
```
- `reference_fa`: path to directory comprising reference genom [*.fa] as well as genome index [*.fai]
- `bmp_manifest`: path to GSA array specific BMP manifest file
- `csv_manifest`: path to GSA array specific CSV manifest file
- `cluster_file`: path to GSA array specific cluster file
- `gsa_idats_dir`: directory comprising GSA idats
- `sample_sheet`: sample sheet contianing at least info about `sample name` and `array position`
- `array_position`: sample sheet column name containing information about array position (sentrix ID and sentrix position) e.g. `205723740073_R03C02`
- `sample_name`: sample sheet column name containing sample name
- `results_dir`: results directory


**IMPORTANT NOTES:** 
- Only IDATs common between `gsa_idats_dir/` and  `sample_sheet` are going to be analysed!
- BMP manifest, CSV manifest as well as cluster file are product specific files therefore should be downloaded from Illumina Product Files page e.g. [GSA v3.0](https://emea.support.illumina.com/downloads/infinium-global-screening-array-v3-0-product-files.html)
- Reference genome along with indexes should be downloaded from [Illumina](https://emea.support.illumina.com/downloads/genome-fasta-files.html)

Workflow config example:

```
mQTL 1.0v Worklow one [SNPs arrays]
==============

Config:
==============
Ref. genome directory [should include *.fa and *.fai files] [--reference_fa <path>]: ../test/resources/
BMP manifest [--bmp_manifest <path>]: ../test/resources/GSA-24v3-0_A1.bpm
CSV manifest [--csv_manifest <path>]: ../test/resources/GSA-24v3-0_A1.csv
Cluster file [--cluster_file <path>]: ../test/resources/GSA-24v3-0_A1_ClusterFile.egt
Number of CPUs [--CPUs <int>]: 10 [default: 10]

PLINK params:
==============
MAF [--MAF <float>]: 0.01 [default: 0.01]
HWE [--HWE <float>]: 1e-50 [default: 1e-50]
GENO [--GENO <float>]: 0.1 [default: 0.1]

Input:
==============
GSA idats [--gsa_idats_dir <path>]: ../test/GSA/
Sample sheet [--sample_sheet <path>]: ../test/sample_sheet.csv
Array position column [--array_position <str>]: ArrayPicker
Sample name column [--sample_name <str>]: Sample_Name

Output:
==============
results directory [--results_dir <path>]: ../test/results
```

Workflow output is placed in <results_dir/flow_one> directory, and includes:
- `filtered_merged.bcf` - merged VCF files, filtered for biallelic loci
- `genotype_table.traw` - BCF file converted to PLINK traw format
- `genotype_table.parquet` - final tabular file comprising information about genotype per sample

### Workflow two: Methylation Array Data Preprocessing

This workflow will preprocess methylation array data from Illumina. It involves the following steps:

```
0. Sanity check
1. Sesame preprocessing: use the Sesame R package for preprocessing (prep code: QCDPB) the methylation data.
3. Cell fraction correction: use linear models to adjust for cellular composition in the samples. [OPTIONAL]
4. Exporting normalized or normalized AND adjusted for tissue composition beta-matrix frame (mynorm.parquet).
5. If cell fraction correction has been performed workflow will also export estimated cell proportions before (CF.csv) and after (CFc.csv) data adjustment.

```

Workflow basic command:

```
nextflow run flow_two.nf --methylation_idats_dir <path> --sample_sheet <path> --array_position <str> --sample_name <str> --results_dir <path>
```

- `methylation_idats_dir`: directory comprising EPIC2/EPIC/450K idats to analyse
- `sample_sheet`: sample sheet contianing at least info about `sample name` and `array position`
- `array_position`: sample sheet column name containing information about array position (sentrix ID and sentrix position) e.g. `205723740073_R03C02`
- `sample_name`: sample sheet column name containing sample name
- `results_dir`: results directory [the same as in flow_one]

Workflow config example:

```
mQTL 1.0v Worklow two [Methylation arrays]
==============

Config:
==============
Cell fraction correction [--correction <boolean: true/false>]: true [default: false]
Deconvolution method [--deconvolution_method <str: CP/RPC/CBS>]: RPC [default: RPC]
Collapse methylation readings to the cg prefix [--collapse_prefix <boolean: true/false>]: true [default: true]
Number of CPUs [--CPUs <int>]: 10 [default: 10]

Input:
==============
Methylation idats [--methylation_idats_dir <path>]: ../test/EPIC/
Sample sheet [--sample_sheet <path>]: ../test/sample_sheet.csv
Array position column [--array_position <str>]: Sentrix_Info
Sample name column [--sample_name <str>]: Sample_Name

Output:
==============
results directory [--results_dir <path>]: ../test/results/
```

Workflow output is placed in <results_dir/flow_two> directory, and includes:
- `mynorm.parquet` - normalized beta-matrix
- `mynorm_corrected.parquet` - normalized beta-matrix corrected for tissue composition (applicable only for blood samples, only if --correction true)
- `CF.csv` - estimated cellular fractions before correction (only if --correction true)
- `CFc.csv`- estimated cellular fractions after correction (only if --correction true)


### Workflow three: mQTL Identification
This workflow will analyze data from flow_one and flow_two to identify mQTLs. It involves the following steps:

```
1. mQTL analysis: Integrate genetic and methylation data to identify mQTLs.
2. Filtering: extracting significant mQTLs based on provided criteria.
3. Annotating: mQTL annotations using VEP webserver.
4. Clumping: Perform clumping to group SNPs based on linkage disequilibrium.
```

Workflow basic command:

```
nextflow run flow_three.nf --csv_gsa_manifest <path> --csv_methylation_manifest <path> --results_dir <path>
```

- csv_gsa_manifest: GSA manifest file
- csv_methylation_manifest: 450K/EPIC/EPIV2 manifest file 
- results_dir: results directory [the same as in flow_one and flow_two]

Workflow config example:
```
mQTL 1.0v Worklow three [mQTL]
==============

Config:
==============
GSA manifest [--csv_gsa_manifest <path>]: ../test/resources/GSA-24v3-0_A1.csv
Methylation manifest [--csv_methylation_manifest <path>]: ../test/resources/infinium-methylationepic-v-1-0-b5-manifest-file.csv
Genome assembly [--genome_assembly <str>]: GRCh37 [default: GRCh37 (important if VEP is ON)]
VEP annotations [--vep_annotations <boolean: true/false>]: false [default: false (very slow)]
Alpha [--alpha <float>]: 0.05 [default: 0.05]
Slope [--slope <float>]: 0.05 [default: 0.05]
Distance [--distance <int>]: 50000 [default: 50000]
Number of CPUs [--CPUs <int>]: 20 [default: 10]

PLINK config:
==============
P1: [--clump_p1 <float>]: 0.0001 [default: 0.0001]
P2: [--clump_p2 <float>]: 0.01 [default: 0.01]
R1: [--clump_r2 <float>]: 0.50 [default: 0.5]
KB: [--clump_kb <int>]: 250 [default: 250]

Input:
==============
Results from flow_one and flow_two are expected to be in [--results_dir <path>]: ../test/results/

Output:
==============
results directory [--results_dir <path>]: ../test/results/
```

Workflow output is placed in <results_dir/flow_three> directory, and includes:
- 'mQTL.parquet' - tabular file comprising mQTL stats
- 'filtered_mQTL.parquet' - tabular file comprising filtered mQTLs based on --alpha and --slope parameters 
- 'cpg_list.txt' and 'rs_list.txt' - lists of CpGs as well as rsIDs for all significant mQTLs
- 'cpg_pval.txt' and 'rs_pval.txt' - lists of all assesed CpGs/SNPs along with FDR corrected p-value (Benjamini/Yekutieli)
- 'report' and 'report_summary.html' - VEP output generated for 'rs_list.txt'
- 'plink.clumped' - clumping results generated based on 'rs_pval.txt' and 'filtered_merged.bcf'
  
### Contributing
Contributions to the project are welcome. Please fork the repository, make your changes, and submit a pull request.

### License
This project is licensed under the MIT License.

### Contact
For questions or issues, please contact [Jan Bińkowski] at [jan.binkowski@pum.edu.pl].
