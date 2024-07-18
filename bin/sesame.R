#!/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Two arguments <idats dir path> and <n CPUs> must be supplied")
} else {
  data_path = args[1]
  n_cpus = args[2]
}

library(BiocParallel)
library(sesame)
library(arrow)

message("Running: ", data_path)
message("n CPUs: ", n_cpus)
message("Loading idats")

betas <- openSesame(data_path, prep="QCDPB", BPPARAM = MulticoreParam(n_cpus))
betas = betas %>% as.data.frame
betas$CpG <- rownames(betas)

betas <- betas[grep("cg", rownames(betas)),]
write_parquet(betas, "mynorm.parquet")

message("Saving file")
message("N CpGs --> ", length(betas$CpG))
message("N samples --> ", length(colnames(betas)) - 1)
message("DONE")
