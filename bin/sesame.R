#!/usr/local/bin/R

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied")
} else {
  data_path = args[1]
  n_cpus = args[2]
}

library(sesame)
library(arrow)
library(retry)

message("Running: ", data_path)
message("n CPUs: ", n_cpus)
message("Loading idats")

betas <- openSesame(data_path, prep="QCDPB", BPPARAM = BiocParallel::MulticoreParam(n_cpus))
betas = betas %>% as.data.frame
betas$CpG <- rownames(betas)

betas <- betas[grep("cg", rownames(betas)),]
write_parquet(betas, file.path(data_path, "mynorm.parquet"))

message("Saving file")
message("N CpGs --> ", length(betas$CpG))
message("N samples --> ", length(colnames(betas)) - 1)
message("DONE")
