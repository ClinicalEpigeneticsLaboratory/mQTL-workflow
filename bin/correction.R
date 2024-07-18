#!/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Two arguments <mynorm path> and <deconvolution method> must be supplied")
} else {
    mynorm_input = args[1]
    deconvolution_method = args[2]
}

library(EpiDISH)
library(arrow)
library(dplyr)

# Function implementation
modified_refBase <- function(beta, method=method){
  data(centDHSbloodDMC.m)
  
  frac.m <- epidish(beta, ref.m = centDHSbloodDMC.m, method = deconvolution_method)
  cellFrac <- frac.m$estF
  
  print(colMeans(cellFrac))
  message(names(which.min(colMeans(cellFrac))), "has smallest cell proportion, all other cell proportions will be corrected by linear regression method.")
  
  lm.o <- lm(t(beta) ~ cellFrac[, -1 * which.min(colMeans(cellFrac))])
  tmp.m <- t(lm.o$res) + rowMeans(beta)

  tmp.m[tmp.m <= 0] <- min(tmp.m[which(tmp.m > 0)])
  tmp.m[tmp.m >= 1] <- max(tmp.m[which(tmp.m < 1)])
  
  frac.m2 <- epidish(tmp.m, ref.m = centDHSbloodDMC.m, method = deconvolution_method)
  cellFrac2 <- frac.m2$estF
  
  return(list(CorrectedBeta = tmp.m, 
              CellFractionBeforeCorrection = cellFrac,
              CellFractionAfterCorrection = cellFrac2))
}

mynorm <- read_parquet(mynorm_input) %>% as.data.frame()
rownames(mynorm) <- mynorm$CpG
mynorm$CpG <- NULL

mynorm <- na.omit(mynorm)
print(head(mynorm))
print(dim(mynorm))
results <- modified_refBase(mynorm)

corrected_mynorm <- as.data.frame(results$CorrectedBeta)
corrected_mynorm$CpG <- row.names(corrected_mynorm)
message("Exporting mynorm: ", dim(corrected_mynorm))

write_parquet(corrected_mynorm, "mynorm_corrected.parquet")
write.csv(results$CellFractionBeforeCorrection, "CF.csv")
write.csv(results$CellFractionAfterCorrection, "CFc.csv")
