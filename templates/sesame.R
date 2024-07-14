args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied")
} else {
  data_path = args[1]  
}

library(sesame)
library(arrow)
library(retry)
n <- 20

message("Running: ", data_path)
message("Loading idats")

sdfs <- retry(openSesame(file.path(data_path, "idats/"), func=NULL, prep="QCDPB", BPPARAM = BiocParallel::MulticoreParam(n)), when="BiocParallel errors", max_tries=10)

if (is.data.frame(sdfs)){
    name <- names(searchIDATprefixes(file.path(data_path, "idats/")))
    
    temp_list <- list()
    temp_list[[name]] <- sdfs %>% as.data.frame
    sdfs <- temp_list
}

message("Converting to betas")
betas = do.call(cbind, BiocParallel::bplapply(sdfs, getBetas, BPPARAM = BiocParallel::MulticoreParam(n)))
betas = betas %>% as.data.frame

betas$CpG <- rownames(betas)

message("Extracting total intensities")
intensities = do.call(cbind, BiocParallel::bplapply(sdfs, totalIntensities, BPPARAM = BiocParallel::MulticoreParam(n)))
intensities = intensities %>% as.data.frame

intensities$CpG <- rownames(intensities)

message("Saving files")
message("N CpGs --> ", length(intensities$CpG))
message("N samples --> ", length(colnames(intensities)) - 1)

betas <- betas[grep("cg", rownames(betas)),]
write_parquet(betas, file.path(data_path, "mynorm.parquet"))

intensities <- intensities[grep("cg", rownames(intensities)),]
write_parquet(intensities, file.path(data_path, "intensity.parquet"))
message("DONE")
