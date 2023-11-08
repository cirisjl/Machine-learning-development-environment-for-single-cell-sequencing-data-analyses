library(ggplot2)
library(SingleR)
library(dplyr)
library(celldex)
library(RColorBrewer)
library("here")
source(here::here('tools/formating/formating.R'))
# source("../../formating/formating.R")


RunSeuratQC <- function(input, output, save_anndata=TRUE, assay='RNA', nFeature_min=200, nFeature_max=0, percent_mt_max=5, percent_rb_min=0, path_of_scrublet_calls=here::here('api/tools/qc/scrublet_calls.tsv'), dims=1:10, regress_cell_cycle=FALSE) {
    srat <- tryCatch(
        LoadSeurat(input),
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
    adata_path <- NULL

    if (!is.null(srat)){
        # If assay if provided by the user, then set default_assy to assay.
        if(assay!='RNA') DefaultAssay(srat) <- assay
        default_assay <- DefaultAssay(srat)
        print(default_assay)

        # Check either the default assay of the Seurat object is "RNA" or the assay is provided by the user.
        if(default_assay=='RNA' | default_assay==assay){
            DefaultAssay(srat) <- assay
            if(!paste0("nCount_", default_assay) %in% names(x = srat[[]])) srat[[paste0("nCount_", default_assay)]] <- colSums(x = srat[[default_assay]], slot = "counts")  # nCount of the default assay
            if(!paste0("nFeature_", default_assay) %in% names(x = srat[[]])) srat[[paste0("nFeature_", default_assay)]] <- colSums(x = GetAssayData(object = srat[[default_assay]], slot = "counts") > 0)  # nFeature of the default assay
            
            # Calculate the percentage of mitocondrial per cell and add to the metadata.
            if(! "percent.mt" %in% names(x = srat[[]])) srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^MT-")
            # Calculate the proportion gene expression that comes from ribosomal proteins.
            if(! "percent.rb" %in% names(x = srat[[]])) srat[["percent.rb"]] <- PercentageFeatureSet(srat, pattern = "^RP[SL]")

            # Percentage hemoglobin genes - includes all genes starting with HB except HBP.
            if(! "percent.hb" %in% names(x = srat[[]])) srat[["percent.hb"]] <- PercentageFeatureSet(srat, pattern = "^HB[^(P)]")
            if(! "percent.plat" %in% names(x = srat[[]])) srat[["percent.plat"]] <- PercentageFeatureSet(srat, pattern = "PECAM1|PF4")

            # Add the doublet annotation
            doublets <- read.table(path_of_scrublet_calls, header = F, row.names = 1)
            colnames(doublets) <- c("Doublet_score", "Is_doublet")
            if(! "Is_doublet" %in% names(x = srat[[]])) {
                srat <- AddMetaData(srat, doublets)
            }

            # print(head(srat@meta.data))

            srat <- subset(srat, subset = paste0("nFeature_", default_assay) > nFeature_min & percent.mt < percent_mt_max)
            if(nFeature_max != 0) srat <- subset(srat, subset = paste0("nFeature_", default_assay) < nFeature_max)
            if(percent_rb_min != 0)  srat <- subset(srat, subset = percent.rb > percent_rb_min)
            srat <- subset(srat, subset = Is_doublet != 'True' | is.na(Is_doublet))
            srat <- NormalizeData(srat, normalization.method = "LogNormalize", scale.factor = 10000)
            srat <- FindVariableFeatures(srat, selection.method = "vst")
            srat <- ScaleData(srat, features = rownames(srat))

            # PCA
            # srat <- RunPCA(srat, features = VariableFeatures(srat), ndims.print = 6:10, nfeatures.print = 10)
            srat <- RunPCA(srat, features = VariableFeatures(srat))

            if(regress_cell_cycle){
                RegressCellCycle(srat)
            }

            srat <- FindNeighbors(srat, dims=dims)
            srat <- FindClusters(srat, resolution = 0.5)
            # TSNE
            srat <- RunTSNE(srat, dims=dims)
            # UMAP
            srat <- RunUMAP(srat, dims=dims)

            assay_names <- names(srat@assays)
            metadata <- srat@meta.data
            nCells <- ncol(srat)
            nGenes <- nrow(srat)
            genes <- rownames(srat)
            cells <- Cells(srat)
            HVGsID <- srat[[assay]]@var.features
            if('pca' %in% names(srat@reductions)) pca <- Embeddings(object = srat, reduction = "pca")
            if('tsne' %in% names(srat@reductions)) tsne <- Embeddings(object = srat, reduction = "tsne")
            if('umap' %in% names(srat@reductions)) umap <- Embeddings(object = srat, reduction = "umap")

            srat@meta.data <- .regularise_df(srat@meta.data, drop_single_values=FALSE, drop_na_values=TRUE)

            SaveH5Seurat(srat, filename=output, overwrite=TRUE, verbose=FALSE)
            print("Seurat object is saved successfully.")
            if(save_anndata){
                adata_path <- Convert(output, dest = "h5ad" , overwrite = TRUE)
                print("AnnData object is saved successfully.")
            }
        } else {
            assay_names <- names(srat@assays)
            default_assay <- DefaultAssay(srat)
        }  
    }
    rm(srat)
    gc()
    list(default_assay=default_assay, assay_names=assay_names, metadata=metadata, nCells=nCells, nGenes=nGenes, genes=genes, cells=cells, HVGsID=HVGsID, pca=pca, tsne=tsne, umap=umap, adata_path=adata_path)
}


# RunSeuratQC <- function(srat, output, nFeature_min=200, nFeature_max=0, percent_mt_max=5, percent_rb_min=0, path_of_scrublet_calls=here::here('api/tools/qc/scrublet_calls.tsv'), dims=1:10, regress_cell_cycle=FALSE) {
#     default_assay <- NULL
#     assay_names <- NULL
#     metadata <- NULL
#     HVGsID <- NULL
#     nGenes <- 0
#     nCells <- 0
#     genes <- NULL
#     cells <- NULL
#     pca <- NULL
#     tsne <- NULL
#     umap <- NULL

#     if (!is.null(srat)){
#         default_assay <- DefaultAssay(srat)
#         # print(default_assay)
        
#         if(!paste0("nCount_", default_assay) %in% names(x = srat[[]])) srat[[paste0("nCount_", default_assay)]] <- colSums(x = srat[[default_assay]], slot = "counts")  # nCount of the default assay
#         if(!paste0("nFeature_", default_assay) %in% names(x = srat[[]])) srat[[paste0("nFeature_", default_assay)]] <- colSums(x = GetAssayData(object = srat[[default_assay]], slot = "counts") > 0)  # nFeature of the default assay
        
#         # Calculate the percentage of mitocondrial per cell and add to the metadata.
#         if(! "percent.mt" %in% names(x = srat[[]])) srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^MT-")
#         # Calculate the proportion gene expression that comes from ribosomal proteins.
#         if(! "percent.rb" %in% names(x = srat[[]])) srat[["percent.rb"]] <- PercentageFeatureSet(srat, pattern = "^RP[SL]")

#         # Percentage hemoglobin genes - includes all genes starting with HB except HBP.
#         if(! "percent.hb" %in% names(x = srat[[]])) srat[["percent.hb"]] <- PercentageFeatureSet(srat, pattern = "^HB[^(P)]")
#         if(! "percent.plat" %in% names(x = srat[[]])) srat[["percent.plat"]] <- PercentageFeatureSet(srat, pattern = "PECAM1|PF4")

#         # Add the doublet annotation
#         doublets <- read.table(path_of_scrublet_calls, header = F, row.names = 1)
#         colnames(doublets) <- c("Doublet_score", "Is_doublet")
#         if(! "Is_doublet" %in% names(x = srat[[]])) {
#             srat <- AddMetaData(srat, doublets)
#         }
#         srat[['Is_doublet']] <- !is.na(srat[['Is_doublet']])

#         # print(head(srat@meta.data))

#         srat <- subset(srat, subset = paste0("nFeature_", default_assay) > nFeature_min & percent.mt < percent_mt_max)
#         if(nFeature_max != 0) srat <- subset(srat, subset = paste0("nFeature_", default_assay) < nFeature_max)
#         if(percent_rb_min != 0)  srat <- subset(srat, subset = percent.rb > percent_rb_min)
#         srat <- subset(srat, subset = Is_doublet != 'True')
#         srat <- NormalizeData(srat, normalization.method = "LogNormalize", scale.factor = 10000)
#         srat <- FindVariableFeatures(srat, selection.method = "vst")
#         srat <- ScaleData(srat, features = rownames(srat))

#         # PCA
#         srat <- RunPCA(srat, features = VariableFeatures(srat), ndims.print = 6:10, nfeatures.print = 10)

#         if(regress_cell_cycle){
#             RegressCellCycle(srat)
#         }

#         srat <- FindNeighbors(srat, dims=dims)
#         srat <- FindClusters(srat, resolution = 0.5)
#         # TSNE
#         srat <- RunTSNE(srat, dims=dims)
#         # UMAP
#         srat <- RunUMAP(srat, dims=dims)

#         assay_names <- names(srat@assays)
#         default_assay <- DefaultAssay(srat)
#         metadata <- srat@meta.data
#         nCells <- ncol(srat)
#         nGenes <- nrow(srat)
#         genes <- rownames(srat)
#         cells <- Cells(srat)
#         HVGsID <- srat[[assay]]@var.features
#         if('pca' %in% names(srat@reductions)) pca <- Embeddings(object = srat, reduction = "pca")
#         if('tsne' %in% names(srat@reductions)) tsne <- Embeddings(object = srat, reduction = "tsne")
#         if('umap' %in% names(srat@reductions)) umap <- Embeddings(object = srat, reduction = "umap")
#     }

#     SaveH5Seurat(srat, filename=output, overwrite=TRUE, verbose=FALSE)
#     print("Seurat object is saved successfully.")
#     rm(srat)
#     gc()
#     list(default_assay=default_assay, assay_names=assay_names, metadata=metadata, nCells=nCells, nGenes=nGenes, genes=genes, cells=cells, HVGsID=HVGsID, pca=pca, tsne=tsne, umap=umap)
# }


RegressCellCycle <- function(srat){
    # Read in the expression matrix The first row is a header row, the first column is rownames
    # exp.mat <- read.table(file = "../tools/qc/nestorawa_forcellcycle_expressionMatrix.txt", header = TRUE,
    #     as.is = TRUE, row.names = 1)

    # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
    # segregate this list into markers of G2/M phase and markers of S phase
    # s.genes <- cc.genes$s.genes
    # g2m.genes <- cc.genes$g2m.genes

    # srat <- CellCycleScoring(srat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
    srat <- CellCycleScoring(srat, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes, set.ident = TRUE)
    srat <- ScaleData(srat, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(srat))

    srat
}