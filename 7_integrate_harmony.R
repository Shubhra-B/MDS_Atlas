library(Seurat)

merged_object_MDS_harmony <- readRDS("merged_object_MDS_renamed_cleaned.rds") |>
  IntegrateLayers(
    method = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction = "harmony",
    verbose = FALSE
  )

merged_object_MDS_harmony[["RNA"]] <- JoinLayers(merged_object_MDS_harmony[["RNA"]])

merged_object_MDS_harmony <- merged_object_MDS_harmony |>
  FindNeighbors(dims = 1:20, reduction = "harmony") |>
  FindClusters(resolution = 2, cluster.name = "harmony_clusters") |>
  RunUMAP(dims = 1:20, reduction = "harmony", reduction.name = "umap.harmony")

saveRDS(merged_object_MDS_harmony, "merged_object_MDS_harmony.rds")

p5 <- DimPlot(
  merged_object_MDS_harmony,
  reduction = "umap.harmony",
  group.by = "sample"
)
print(p5)
dev.off()
