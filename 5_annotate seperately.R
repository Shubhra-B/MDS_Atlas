#I like to look at the samples first seperately and then merge,integrate


# Load reference once (avoid reloading in each loop)
ref <- loadHDF5SummarizedExperiment(
  dir = "NovershternHematopoieticData",
  prefix = ""
)

for (file in list.files(pattern = "_seur_filt\\.rds$")) {
  base_name <- sub("\\_seur_filt$", "", file)  
  obj <- readRDS(file)
  obj <- NormalizeData(obj, normalization.method = "LogNormalize")
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  Identify the 10 most highly variable genes
  top10 <- head(VariableFeatures(obj), 10)
  all.genes <- rownames(obj)
  obj <- ScaleData(obj, features = all.genes, verbose = FALSE)	

  obj <- RunPCA(obj, features = VariableFeatures(object = obj), verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:10)
  obj <- FindClusters(obj, resolution = 0.5)
  obj <- RunUMAP(obj, dims = 1:10)
  
  pdf(paste0(base_name, "_UMAP.pdf"), width = 6, height = 5)
  print(DimPlot(obj, label = TRUE))
  dev.off()
  
  # SingleR annotation
  pred <- SingleR(
    test = GetAssayData(obj, slot = "data"),
    ref = ref,
    labels = ref$label.main
  )
  obj@meta.data$predictions_nov <- pred$labels
  saveRDS(obj, paste0(base_name, "_proc.rds"))
  
  cat("Saved processed object:", paste0(base_name, "_annotated.rds"), "\n")
}
