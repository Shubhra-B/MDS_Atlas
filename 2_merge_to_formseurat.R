library(Seurat)

#all the files have barcodes,mtx and features all seperately we need to have seurat object
for (dir in list.dirs(".", recursive = FALSE)) {
  
  data <- Read10X(
    data.dir = dir,
    gene.column = 2,
    cell.column = 1)
  
  # Create Seurat object
  srat <- CreateSeuratObject(data, project = base_name, assay = "RNA")
 saveRDS(srat, paste0(base_name, "_seur.rds")) 
  cat("Saved", paste0(base_name, "_seur.rds"), "\n")}
