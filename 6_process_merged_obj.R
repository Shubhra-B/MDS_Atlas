#start with the merged object here, I had to do some manual shit to work with the label
library(Seurat)

merged_object_MDS_renamed_cleaned <- readRDS("merged_object_MDS_renamed_cleaned.rds") |>
  NormalizeData() |>
  FindVariableFeatures() |>  # initial run
  FindVariableFeatures(nfeatures = 5000) |>  # override with 5000 features
  ScaleData(features = rownames(.)) |>
  RunPCA() |>
  FindNeighbors(dims = 1:20) |>
  FindClusters() |>
  RunUMAP(dims = 1:20)
  ElbowPlot(merged_object_MDS_renamed_cleaned)
  p1 <- DimPlot(merged_object_MDS_renamed_cleaned, group.by = "sample")
  print(p1)
  dev.off()

saveRDS(merged_object_MDS_renamed_cleaned, "merged_object_MDS_renamed_cleaned.rds")
