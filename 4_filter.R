library(Seurat)
library(ggplot2)

# Loop over only *_seur.rds files
for (file in list.files(pattern = "_seur\\.rds$")) {
  
  base_name <- sub("\\.rds$", "", file)   # remove .rds extension
  
  x <- readRDS(file)
  # I use this quite loose but we go easy first and readjust later
  x <- subset(x, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent_mito < 10)

  # Save filtered seur object
  saveRDS(x,paste0(base_name, "_filt.rds"))

  cat("filtered the object", base_name, "\n")
}
