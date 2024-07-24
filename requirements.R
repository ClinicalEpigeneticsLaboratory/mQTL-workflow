install.packages(c("BiocManager", "dplyr", "arrow"))
BiocManager::install(c("sesame", "EpiDISH"))

library(sesame)
sesameDataCache()
