# Function to perform cell-fraction correction based on ChAMP implentation and RPC method from
# EpiDISH package.

library("EpiDISH")
library("arrow")
library("dplyr")

args = commandArgs(trailingOnly=TRUE)
message("Input: ", args[1], "\n",
        "output mynorm: ", args[2], "\n", 
        "output cell proprotions estimation: ", args[3],
        "\n")

mynorm_input = args[1]
mynorm_output = args[2]
cfp_output_dir = args[3]

# Function implementation
modified_refBase <- function(beta, method="RPC"){
  
  data(centDHSbloodDMC.m)
  
  frac.m <- epidish(beta, ref.m = centDHSbloodDMC.m, method = method)
  cellFrac <- frac.m$estF
  
  print(colMeans(cellFrac))
  message(names(which.min(colMeans(cellFrac))), "has smallest cell proportion, all other cell proportions will be corrected by linear regression method.")
  
  lm.o <- lm(t(beta) ~ cellFrac[, -1 * which.min(colMeans(cellFrac))])
  tmp.m <- t(lm.o$res) + rowMeans(beta)

  tmp.m[tmp.m <= 0] <- min(tmp.m[which(tmp.m > 0)])
  tmp.m[tmp.m >= 1] <- max(tmp.m[which(tmp.m < 1)])
  
  frac.m2 <- epidish(tmp.m, ref.m = centDHSbloodDMC.m, method = method)
  cellFrac2 <- frac.m2$estF
  
  return(list(CorrectedBeta = tmp.m, 
              CellFractionBeforeCorrection = cellFrac,
              CellFractionAfterCorrection = cellFrac2))
}

mynorm <- read_parquet(mynorm_input) %>% as.data.frame()
message("Loaded mynorm: ", dim(mynorm))

rownames(mynorm) <- mynorm$CpG
mynorm$CpG <- NULL

results <- modified_refBase(mynorm)

corrected_mynorm <- as.data.frame(results$CorrectedBeta)
corrected_mynorm$CpG <- row.names(corrected_mynorm)
message("Exporting mynorm: ", dim(corrected_mynorm))

write_parquet(corrected_mynorm, mynorm_output)
write.csv(results$CellFractionBeforeCorrection, file.path(cfp_output_dir, "CF.csv"))
write.csv(results$CellFractionAfterCorrection, file.path(cfp_output_dir, "CFc.csv"))
