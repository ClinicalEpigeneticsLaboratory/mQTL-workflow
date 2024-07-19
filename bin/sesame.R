#!/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Two arguments <idats dir path> and <n CPUs> must be supplied")
} else {
  data_path = args[1]
  n_cpus = args[2]
  collapse_prefix = args[3]
}

library(BiocParallel)
library(sesame)
library(arrow)

message("Running: ", data_path)
message("n CPUs: ", n_cpus)
message("Collapse prefix: ", collapse_prefix)
message("Loading idats")

if (collapse_prefix == "true"){
    betas <- openSesame(data_path, prep="QCDPB", BPPARAM = MulticoreParam(n_cpus), collapseToPfx = TRUE, collapseMethod = "mean")
} else {
    betas <- openSesame(data_path, prep="QCDPB", BPPARAM = MulticoreParam(n_cpus), collapseToPfx = FALSE)
}

betas = betas %>% as.data.frame
betas$CpG <- rownames(betas)

betas <- betas[grep("cg", rownames(betas)),]
write_parquet(betas, "mynorm.parquet")

message("Saving file")
message("N CpGs --> ", length(betas$CpG))
message("N samples --> ", length(colnames(betas)) - 1)
message("DONE")
