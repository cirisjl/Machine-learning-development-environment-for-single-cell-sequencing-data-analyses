library(ggplot2)
library(SingleR)
library(dplyr)
library(celldex)
library(RColorBrewer)
source(here::here('api/tools/formating/formating.R'))


RunSeuratQC <- function(input_path, output_path, assay='RNA', regress_cell_cycle=FALSE) {
    srat <- tryCatch(
        LoadSeurat(path),
        error = function(e) {
            stop("The file format is not supported.")
            print(e)
        }
    )

    default_assay <- NULL
    assay_names <- NULL
    metadata <- NULL
    HVGsID <- NULL
    nGenes <- 0
    nCells <- 0
    genes <- NULL
    cells <- NULL
    pca <- NULL
    tsne <- NULL
    umap <- NULL

    if (!is.null(srat)){
        assay_names <- names(srat@assays)
        default_assay <- DefaultAssay(srat)
    }

    if (!is.null(srat) && DefaultAssay(srat)=='RNA'){
        srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^MT-")
        srat <- subset(srat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
        srat <- NormalizeData(srat, normalization.method = "LogNormalize", scale.factor = 10000)
        srat <- FindVariableFeatures(srat, selection.method = "vst")
        srat <- ScaleData(srat, features = rownames(srat))

        srat <- RunPCA(srat, features = VariableFeatures(srat), ndims.print = 6:10, nfeatures.print = 10)

        if(regress_cell_cycle){
            RegressCellCycle(srat)
        }

        metadata <- srat@meta.data
        nCells <- ncol(srat)
        nGenes <- nrow(srat)
        genes <- rownames(srat)
        cells <- Cells(srat)
        HVGsID <- srat[[assay]]@var.features
        if('pca' %in% names(srat@reductions)) pca <- Embeddings(object = srat, reduction = "pca")
        if('tsne' %in% names(srat@reductions)) tsne <- Embeddings(object = srat, reduction = "tsne")
        if('umap' %in% names(srat@reductions)) umap <- Embeddings(object = srat, reduction = "umap")
    }

    SaveH5Seurat(srat, filename=output_path, overwrite=TRUE, verbose=FALSE)
    print("Seurat object is saved successfully.")
    rm(srat)
    gc()
    list(default_assay=default_assay, assay_names=assay_names, metadata=metadata, nCells=nCells, nGenes=nGenes, genes=genes, cells=cells, HVGsID=HVGsID, pca=pca, tsne=tsne, umap=umap)
}


RegressCellCycle <- function(srat){
    # Read in the expression matrix The first row is a header row, the first column is rownames
    # exp.mat <- read.table(file = "../tools/qc/nestorawa_forcellcycle_expressionMatrix.txt", header = TRUE,
    #     as.is = TRUE, row.names = 1)

    # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
    # segregate this list into markers of G2/M phase and markers of S phase
    # s.genes <- cc.genes$s.genes
    # g2m.genes <- cc.genes$g2m.genes

    # srat <- CellCycleScoring(srat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
    srat <- CellCycleScoring(srat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
    srat <- ScaleData(srat, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(srat))

    srat
}