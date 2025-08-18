
library(dplyr)
library(pheatmap)
library(ggplot2)
library(tibble)
library(tidyr)
library(Seurat)
library(hdf5r)
library(sctransform)
library(SummarizedExperiment)
library(patchwork)
library(celldex)
library(SingleR)
library(HDF5Array)
library("BiocFileCache")
library(matrixStats)

packageVersion("Seurat")
packageVersion("SeuratObject")
library("BiocFileCache")
library(matrixStats)


cleaned_object_seurat<-readRDS("cleaned_object_seurat.rds")

#make sure you have the NovershternHematopoieticData installed by https://bioconductor.org/packages/release/data/experiment/vignettes/celldex/inst/doc/userguide.html 
ref=loadHDF5SummarizedExperiment(dir="NovershternHematopoieticData", prefix="") #need to have BiocFileCache at 2.11.1  or will throw error #https://stat.ethz.ch/pipermail/bioc-devel/2023-October/020003.html/pipermail/bioc-devel/2023-October/020003.html


#need to have matrixStats 1.1.0 or will throw error 
predictions_nov <- SingleR(test=GetAssayData(cleaned_object_seurat,layers = "data"),ref=ref, labels=ref$label.main)
#re-assign the label
cleaned_object_seurat@meta.data$predictions_nov <- predictions_nov$labels
saveRDS(cleaned_object_seurat,"cleaned_object_seurat_nov_annotated.rds")

#look at the umap
dp3=DimPlot(cleaned_object_seurat_nov_annotated,reduction="umap.harmony",group.by="predictions_nov")
dp3
dev.off()
