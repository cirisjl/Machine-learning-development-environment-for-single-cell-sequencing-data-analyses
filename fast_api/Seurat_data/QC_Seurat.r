library(Seurat)
library(ggplot2)
library(SingleR)
library(dplyr)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)
library(SeuratDisk)
library(SeuratData)
library(anndata)
library(patchwork)
src <- "D:/Master's/ML_backend_AISingleCell/Machine-learning-development-environment-for-single-cell-sequencing-data-analyses/fast_api/Seurat_data/formatting.R"
source(src)

input_path <- "D:/Master's/ML_backend_AISingleCell/Machine-learning-development-environment-for-single-cell-sequencing-data-analyses/fast_api/Seurat_data/input_data"
output_path <- "D:/Master's/ML_backend_AISingleCell/Machine-learning-development-environment-for-single-cell-sequencing-data-analyses/fast_api/Seurat_data/output_seurat/output_seurat"
input_filename <- paste(input_path, "/scrublet_calls.tsv", sep="")

    
input <- input_path
output_seurat <- output_path
srat <- load_seurat(input)
str(srat)

meta <- srat@meta.data
dim(meta)

head(meta)

summary(meta$nCount_RNA)

summary(meta$nFeature_RNA)

srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^MT-")

srat[["percent.rb"]] <- PercentageFeatureSet(srat, pattern = "^RP[SL]")


doublets <- read.table(input_filename,header = F,row.names = 1)
colnames(doublets) <- c("Doublet_score","Is_doublet")
srat <- AddMetaData(srat,doublets)
head(srat[[]])

VlnPlot(srat, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"),ncol = 4,pt.size = 0.1) &
theme(plot.title = element_text(size=10))

FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "percent.mt")

FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "percent.rb")

FeatureScatter(srat, feature1 = "percent.rb", feature2 = "percent.mt")

FeatureScatter(srat, feature1 = "nFeature_RNA",
               feature2 = "Doublet_score")

srat[['QC']] <- ifelse(srat@meta.data$Is_doublet == 'True','Doublet','Pass')
srat[['QC']] <- ifelse(srat@meta.data$nFeature_RNA < 500 & srat@meta.data$QC == 'Pass','Low_nFeature',srat@meta.data$QC)
srat[['QC']] <- ifelse(srat@meta.data$nFeature_RNA < 500 & srat@meta.data$QC != 'Pass' & srat@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',srat@meta.data$QC,sep = ','),srat@meta.data$QC)
srat[['QC']] <- ifelse(srat@meta.data$percent.mt > 15 & srat@meta.data$QC == 'Pass','High_MT',srat@meta.data$QC)
srat[['QC']] <- ifelse(srat@meta.data$nFeature_RNA < 500 & srat@meta.data$QC != 'Pass' & srat@meta.data$QC != 'High_MT',paste('High_MT',srat@meta.data$QC,sep = ','),srat@meta.data$QC)
table(srat[['QC']])


VlnPlot(subset(srat, subset = QC == 'Pass'),
		features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 4, pt.size = 0.1) &
theme(plot.title = element_text(size=10))


SaveH5Seurat(srat, filename = output_seurat, overwrite = TRUE)
# Convert(srat, dest = "h5ad", overwrite = TRUE)