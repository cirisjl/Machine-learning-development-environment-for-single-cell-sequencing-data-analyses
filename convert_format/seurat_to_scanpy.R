library(qs)
library(Seurat)
library(SeuratDisk)

setwd('C:/Users/flyku/Documents/GitHub/Machine-learning-development-environment-for-single-cell-sequencing-data-analyses/benchmark')
set.seed(42)

# Read qs data
combined <- qs::qread('./test.qsave')

# Save qs data
qs::qsave(combined, './test.qsave')

# Set default assay
Idents(combined) <- combined$seurat_clusters
DefaultAssay(combined) <- "RNA"

# Export Seurat to AnnData
SaveH5Seurat(combined, filename = paste0("test_raw.h5Seurat"),  overwrite = T)
Convert(paste0("test_raw.h5Seurat"), dest = "h5ad", overwrite = T)

# Import AnnData to Seurat
Convert("test.h5ad", dest = "h5seurat", overwrite = TRUE)
test_object <- LoadH5Seurat("test.h5seurat")
test_object

