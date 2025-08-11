library(Seurat)
library(ggplot2)

for (file in list.files(pattern = "_seur\\.rds$")) {
  base_name <- sub("\\.rds$", "", file)   # remove .rds extension
  x <- readRDS(file)  
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
    p <- VlnPlot(x,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3)
  pdf(paste0(base_name, "_QC.pdf"), width = 10, height = 4)
  print(p)
  dev.off()
  
  cat("Saved QC plot for", base_name, "\n")}
