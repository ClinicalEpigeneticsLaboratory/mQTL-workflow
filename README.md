# mQTL-workflow
This Nextflow project is designed to identify methylation quantitative trait loci (mQTL), which are interactions between the genome and the methylome. The project consists of three main workflows:

1. `flow_one` - Preprocess SNP arrays data (Illumina).
2. `flow_two` - Preprocess methylation arrays data (Illumina).
3. `flow_three` - Identify mQTLs.

### Workflows Overview

#### Workflow two: SNP Array Data Preprocessing
The first workflow preprocesses SNP array data from Illumina microarrays. It involves the following steps:

0. Sanity check (not implemented yet)
1. Convert idat to gtc: Convert raw idat files to gtc format.
2. Convert gtc to vcf: Convert gtc files to VCF format.
3. Index and merge VCF files: Index and merge VCF files into a BCF file.
4. Filtering and formatting:
    1. Use PLINK to filter the merged BCF file.
    2. Apply filters: minor allele frequency (--maf 0.05), Hardy-Weinberg equilibrium (--hwe 1e-50), genotype missingness (--geno 0.1), and autosome restriction.
    3. Recode to a transposed format.
5. Exporting results files (merged.bcf, *.frq as well as *.traw).

Workflow Params:

```
mQTL v1.0 Worklow one [SNPs arrays]
==============

Config:
==============
Reference genome: ~/projects/mQTL/GSA/resources/GRCh37_genome.fa
BMP manifest: ~/projects/mQTL/GSA/resources/GSA-24v3-0_A1.bpm
CSV manifest: ~/projects/mQTL/GSA/resources/GSA-24v3-0_A1.csv
Cluster file: ~/projects/mQTL/GSA/resources/GSA-24v3-0_A1_ClusterFile.egt
Number of CPUs: 20

Input:
==============
GSA idats: GSA/

Output:
==============
results directory: results/
```

#### Workflow two: Methylation Array Data Preprocessing
Note: This workflow is not yet implemented.

This workflow will preprocess methylation array data from Illumina. Planned steps include:

0. Sanity check
1. Sesame preprocessing: use the Sesame R package for preprocessing the methylation data.
2. Cell fraction correction: use linear models to adjust for cellular composition in the samples.

#### Workflow three: mQTL Identification
Note: This workflow is not yet implemented.

This workflow will analyze data from flow_one and flow_two to identify mQTLs. Planned steps include:

1. mQTL analysis: Integrate genetic and methylation data to identify mQTLs.
2. Clumping: Perform clumping to group SNPs based on linkage disequilibrium.
3. Enrichment analysis: Conduct enrichment analysis to interpret the biological significance of the identified mQTLs.


### Getting Started

#### Docker Container
A Docker container will be added to the project to standardize the environment and ensure reproducibility.

#### Prerequisites
Ensure you have the following software installed:

- Nextflow
- Docker

#### Running the Workflows
To run the workflows, use the following Nextflow command:


```

nextflow run <workflow_name> -with-apptainer <workflow-specific-params>

```

#### Workflow Configuration
Adjust the configuration parameters in nextflow.config to suit your data and environment.

### Contributing
Contributions to the project are welcome. Please fork the repository, make your changes, and submit a pull request.

### License
This project is licensed under the MIT License.

### Contact
For questions or issues, please contact [Jan Bińkowski] at [jan.binkowski@pum.edu.pl].
