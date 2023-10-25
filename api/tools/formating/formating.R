library(scater)
library(anndata)
library(Seurat)
library(SingleCellExperiment)
library(SeuratDisk)
library(SeuratData)
library(patchwork)
library(Signac)
# library(loomR)

# Get suffix of the file
GetSuffix <- function(path) {
    filename <- basename(path)
    parts <- strsplit(filename,".",fixed = TRUE)
    nparts <- length(parts[[1]])
    return(parts[[1]][nparts])
}

# Read csv/xlsx/h5ad/hdf5/h5/loom/mtx/txt/tab/data/gz file to create AnnData object
LoadAnndata <- function(path) {
    adata <- NULL
    suffix <- tolower(GetSuffix(path))
    if(suffix == "csv"){
        adata <- read_csv(path)
    } else if(suffix == "xlsx"){
        adata <- read_excel(path)
    } else if(suffix == "h5ad"){
        adata <- read_h5ad(path)
    } else if(suffix == "hdf5" || suffix == "h5"){
        adata <- read_hdf(path)
    } 
    # else if(suffix == "loom"){
    #     ad <- read_loom(path)
    # }
     else if(suffix == "mtx"){
        adata <- read_mtx(path)
    } else if(suffix == "txt" || suffix == "tab" || suffix == "data"){
        delim <- DetectDelim(path)
        adata <- read_text(path, delimiter = delim)
    } else if(suffix == "gz"){
        adata <- read_umi_tools(path)
    } else if(suffix == "h5Seurat" || suffix == "h5seurat"){
        Convert(path, dest = "h5ad", overwrite = TRUE, verbose = FALSE)
        adata <- read_h5ad(adata_path)
    } else if(suffix == "rds" || suffix == "robj"){
        srat <- LoadSeurat(path)
        seurat_path <- paste0(tools::file_path_sans_ext(path), ".h5Seurat")
        SaveH5Seurat(srat, filename = seurat_path, overwrite = TRUE, verbose = FALSE)        
        Convert(paste0(tools::file_path_sans_ext(path), ".h5Seurat"), dest = "h5ad" , overwrite = TRUE, verbose = FALSE)
        adata_path <- Convert(seurat_path, dest = "h5ad" , overwrite = TRUE, verbose = FALSE)
        adata <- read_h5ad(adata_path)    } 
    # else if(suffix == "loom"){
    #     srat <- LoadSeurat(path)
    #     SaveH5Seurat(srat, overwrite = TRUE)
    #     Convert(paste0(tools::file_path_sans_ext(path), ".h5Seurat"), dest = "h5ad")
    #     ad <- read_h5ad(paste0(tools::file_path_sans_ext(path), ".h5ad"))
    # } 
    py_to_r_ifneedbe(adata)
}

LoadExpressionMatrix <- function(path) {
    expression_matrix <- NULL
    if(file_test("-d", path)) {
        if(file.exists(file.path(path,"barcodes.tsv")) && file.exists(file.path(path,"genes.tsv")) && file.exists(file.path(path,"matrix.mtx"))){
            expression_matrix <- Read10X(data.dir = path)
        } else if(file.exists(file.path(path,"barcodes.tsv.gz")) && file.exists(file.path(path,"genes.tsv.gz")) && file.exists(file.path(path,"matrix.mtx.gz"))){
            expression_matrix <- Read10X(data.dir = path)
        } else if(file.exists(file.path(path,"count_matrix.mtx.gz")) && file.exists(file.path(path,"features.tsv.gz")) && file.exists(file.path(path,"barcodes.tsv.gz"))){
            expression_matrix <- ReadMtx(mtx = "count_matrix.mtx.gz", features = file.path(path,"features.tsv.gz"), cells = file.path(path,"barcodes.tsv.gz"))
        } else if(file.exists(file.path(path,"count_matrix.mtx")) && file.exists(file.path(path,"features.tsv")) && file.exists(file.path(path,"barcodes.tsv"))){
            expression_matrix <- ReadMtx(mtx = "count_matrix.mtx", features = file.path(path,"features.tsv"), cells = file.path(path,"barcodes.tsv"))           
        } else if(file.exists(file.path(path,"molecules.txt")) && file.exists(file.path(path,"annotation.txt"))){
            delim <- DetectDelim(path)
            molecules <- read.delim(file.path(path,"molecules.txt"), sep = delim, row.names = 1) 
            expression_matrix <- as.matrix(molecules)} 
    } else{       
        suffix <- tolower(GetSuffix(path))
        if(suffix == "h5"){
            expression_matrix  <- Read10X_h5(path, use.names = T)
            if('Gene Expression' %in% names(expression_matrix)) expression_matrix <- expression_matrix$`Gene Expression`
        }
        else if(suffix == "csv"){
            expression_matrix  <- as.matrix(read.csv(path, header=TRUE, row.names=1))
        } 
    }
    expression_matrix
}

LoadSeurat <- function(path, project = NULL) {
    srat <- NULL
    suffix <- tolower(GetSuffix(path))

    if(suffix == "h5Seurat" || suffix == "h5seurat"){
        print("Inside LoadSeurat")
        srat <- LoadH5Seurat(path)
    } else if(suffix == "h5ad"){
        Convert(path, "h5seurat", overwrite = TRUE, assay = "RNA")
        srat <- LoadH5Seurat(paste0(tools::file_path_sans_ext(path), ".h5seurat"))
    } else if(suffix == "rds"){
        robj <- readRDS(path)
        if(class(robj) == 'Seurat'){
            srat <- CreateSeuratObject(counts=robj[['RNA']]@counts, meta.data=robj@meta.data, project = Project(robj))
        } else if(class(robj) == 'SingleCellExperiment'){
            # srat <- as.Seurat(sce, slot = "counts", data = NULL)
            srat <- as.Seurat(robj, slot = "counts")
        }
        rm(robj)
    } else if(suffix == "robj"){
        srat_v2 <- get(load(path))
        if(class(srat_v2) == 'seurat'){
        srat_v3 <- UpdateSeuratObject(srat_v2)
        srat <- CreateSeuratObject(counts=srat_v3[['RNA']]@counts, meta.data=srat_v3@meta.data, project = Project(srat_v3))
        rm(srat_v3)
        }
        rm(srat_v2)
    } else {
        expression_matrix <- LoadExpressionMatrix(path)
        if(!is.null(expression_matrix) && !is.null(project)) {
            srat <- CreateSeuratObject(counts = expression_matrix, project = project)
        } else if (!is.null(expression_matrix)){
            srat <- CreateSeuratObject(counts = expression_matrix)
        }
        rm(expression_matrix) # Erase expression_matrix from memory to save RAM
    }
    # else if(suffix == "loom"){
    #     loom <- connect(filename = path, mode = "r")
    #     srat <- as.Seurat(loom)
    # }
    gc()
    srat
}


LoadSCE <- function(path) {
    sce <- NULL
    if(file_test("-d", path)){
        print("Inside LoadSeurat 3")
        if(file.exists(file.path(path,"molecules.txt")) && file.exists(file.path(path,"annotation.txt"))){
            delim <- DetectDelim(path)
            molecules <- read.delim(file.path(path,"molecules.txt"), sep = delim, row.names = 1) 
            annotation <- read.delim(file.path(path,"annotation.txt"), sep = delim, stringsAsFactors = T)
            sce <- SingleCellExperiment(assays = list(counts = as.matrix(molecules)), colData = annotation)
        } else{
            print("Inside LoadSeurat 4")
            srat <- LoadSeurat(path)
            if(!is.null(srat)){
                sce <- as.SingleCellExperiment(srat)
                rm(srat) # Erase expression_matrix from memory to save RAM
            }
        }
    } else{
        print("Inside LoadSeurat 5")
        suffix <- tolower(GetSuffix(path))
        if(suffix == "rds"){
            sce <- readRDS(path)
            print("Inside LoadSeurat 6")
        } else{
            print("Inside LoadSeurat 7")
            srat <- LoadSeurat(path)
            if(!is.null(srat)){
                sce <- as.SingleCellExperiment(srat)
                rm(srat) # Erase expression_matrix from memory to save RAM
            }
        }
    }
    sce
}


# Return metadata and a list of assays if the default assay is not "RNA"
GetMetadataFromSeurat <- function(path, assay='RNA') {
    srat <- LoadSeurat(path)
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
    # suffix <- tolower(GetSuffix(path))

    if(!is.null(srat)){
        assay_names <- names(srat@assays)
        if(assay != 'RNA' && assay %in% assay_names) DefaultAssay(srat) <- assay
        default_assay <- DefaultAssay(srat)
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
    srat <- NULL
    
    list(default_assay=default_assay, assay_names=assay_names, metadata=metadata, nCells=nCells, nGenes=nGenes, genes=genes, cells=cells, HVGsID=HVGsID, pca=pca, tsne=tsne, umap=umap)
}


# Return a list of assays if the default assay is not "RNA"
ConvertSeuratSCEtoAnndata <- function(path, assay = NULL) {
    assay_names <- NULL
    default_assay <- NULL
    anndata_path <- NULL
    suffix <- tolower(GetSuffix(path))

    srat <- LoadSeurat(path)
    default_assay <- DefaultAssay(srat)
    assay_names <- names(srat@assays)

    if(suffix == "h5Seurat" || suffix == "h5seurat"){
        if(!is.null(assay) && assay %in% assay_names && assay != 'RNA') {
            DefaultAssay(srat) <- assay
            SaveH5Seurat(srat, filename = path, overwrite = TRUE, verbose = FALSE)
            seurat_path <- paste0(tools::file_path_sans_ext(path), "_", assay, ".h5seurat")
            adata_path <- Convert(seurat_path, dest = "h5ad", assay=assay, overwrite = TRUE, verbose = FALSE)
        } else {
            adata_path <- Convert(path, dest = "h5ad", assay=assay, overwrite = TRUE, verbose = FALSE)
        }
    } else if(suffix == "rds" || suffix == "robj"){
        seurat_path <- paste0(tools::file_path_sans_ext(path), ".h5Seurat")
        SaveH5Seurat(srat, filename = seurat_path, overwrite = TRUE, verbose = FALSE)
        adata_path <- Convert(seurat_path, dest = "h5ad" , overwrite = TRUE, verbose = FALSE)
    }
    srat <- NULL

    list(default_assay=default_assay, assay_names=assay_names, anndata_path=anndata_path) # First, check if the anndata_path is NULL; Second, check the assay_names
}


ConvertToAnndata <- function(path, assay = 'RNA') {
    adata_path <- NULL
    suffix <- tolower(GetSuffix(path))
    srat <- LoadSeurat(path)
    seurat_path <- paste0(tools::file_path_sans_ext(path), "_", assay, ".h5seurat")
    if(suffix == "h5Seurat" || suffix == "h5seurat"){
        if(assay != 'RNA') {
            DefaultAssay(srat) <- assay
            SaveH5Seurat(srat, filename = seurat_path, overwrite = TRUE, verbose = FALSE)
            adata_path <- Convert(seurat_path, dest = "h5ad", assay=assay, overwrite = TRUE, verbose = FALSE)
        } else {
            adata_path <- Convert(path, dest = "h5ad", assay=assay, overwrite = TRUE, verbose = FALSE)
        }
    } else if(suffix == "rds"){
        srat <- LoadSeurat(path)
        seurat_path <- paste0(tools::file_path_sans_ext(path), ".h5Seurat")
        SaveH5Seurat(srat, filename = seurat_path, overwrite = TRUE, verbose = FALSE)
        adata_path <- Convert(seurat_path, dest = "h5ad" , overwrite = TRUE, verbose = FALSE)
    } 
    adata_path
}



# save_as_anndata <- function(srat, path, assay = 'RNA'){
#     anndata_path <- gsub("h5seurat", "h5ad", path, ignore.case = TRUE)
#     DefaultAssay(srat) <- assay
#     srat = DietSeurat(
#         srat,
#         counts = TRUE, # so, raw counts save to adata.layers['counts']
#         data = TRUE, # so, log1p counts save to adata.X when scale.data = False, else adata.layers['data']
#         scale.data = FALSE, # if only scaled highly variable gene, the export to h5ad would fail. set to false
#         features = rownames(srat), # export all genes, not just top highly variable genes
#         assays = assay,
#         dimreducs = c("pca","umap"),
#         graphs = c("RNA_nn", "RNA_snn"), # to RNA_nn -> distances, RNA_snn -> connectivities
#         misc = TRUE
#         )
#     if(!MuDataSeurat::WriteH5AD(srat, anndata_path, assay=assay)){
#         anndata_path <- NULL
#     }
#     srat <- NULL
#     anndata_path
# }


#' Automatically detect delimiters in a text file
#'
#' This helper function was written expressly for \code{\link{set_physical}} to
#' be able to automate its \code{recordDelimiter} argument.
#'
#' @param path (character) File to search for a delimiter
#' @param nchar (numeric) Maximum number of characters to read from disk when
#' searching
#'
#' @return (character) If found, the delimiter, it not, \\r\\n
DetectDelim <- function(path, nchar = 1e3) {
  # only look for delimiter if the file exists
  if (file.exists(path)) {
    # readChar() will error on non-character data so
    chars <- tryCatch(
      {
        readChar(path, nchar)
      },
      error = function(e) {
        NA
      }
    )
    search <- regexpr("[,|\\t|;||]+", chars, perl = TRUE)

    if (!is.na(search) && search >= 0) {
      return(substr(chars, search, search + attr(search, "match.length") - 1))
    }
  }
  # readChar() will error on non-character data 
}

py_to_r_ifneedbe <- function(x) {
    if (inherits(x, "python.builtin.object")) {
        py_to_r(x)
    } else {
        x
    }
}

SeuratToCSV <- function(srat, srat_path, assay = 'RNA', slot = "counts"){
    if(assay != 'RNA') slot ="data"
    csv_path <- gsub(".h5Seurat", paste("_", assay, ".csv", sep = ""), srat_path)
    if(!file.exists(csv_path)){
        write.table(as.matrix(GetAssayData(object = srat, assay = assay, slot = slot)), 
        csv_path, sep = ',', row.names = T, col.names = T, quote = F)
    } else{
        print("CSV file already exists.")
    }
    csv_path
}


PlotIntegratedClusters <- function (srat) { 
  ## take an integrated Seurat object, plot distributions over orig.ident
  library(Seurat)
  library(patchwork)
  library(ggplot2)
  library(reshape2)
  library(RColorBrewer)
  
  
  count_table <- table(srat@meta.data$seurat_clusters, srat@meta.data$orig.ident)
  count_mtx   <- as.data.frame.matrix(count_table)
  count_mtx$cluster <- rownames(count_mtx)
  melt_mtx    <- melt(count_mtx)
  melt_mtx$cluster <- as.factor(melt_mtx$cluster)

  cluster_size   <- aggregate(value ~ cluster, data = melt_mtx, FUN = sum)
  
  sorted_labels <- paste(sort(as.integer(levels(cluster_size$cluster)),decreasing = T))
  cluster_size$cluster <- factor(cluster_size$cluster,levels = sorted_labels)
  melt_mtx$cluster <- factor(melt_mtx$cluster,levels = sorted_labels)
  
  colnames(melt_mtx)[2] <- "dataset"
  
  
  p1 <- ggplot(cluster_size, aes(y= cluster,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
    theme_bw() + scale_x_log10() + xlab("Cells per cluster, log10 scale") + ylab("")
  p2 <- ggplot(melt_mtx,aes(x=cluster,y=value,fill=dataset)) + 
    geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
    scale_fill_brewer(palette = "Set2") +
    ylab("Fraction of cells in each dataset") + xlab("Cluster number") + theme(legend.position="top")
  
  p2 + p1 + plot_layout(widths = c(3,1))
  }


LoadMetadata <- function(seurat_obj) {
    metadata <- list()  # Create an empty list to hold metadata
    
    # Get the Default Assay
    metadata$default_assay <- DefaultAssay(object = seurat_obj)
    
    # Get the names of all assays
    metadata$assay_names <- names(seurat_obj@assays)
    if (is.null(metadata$assay_names)) {
        # Return an empty array (vector)
        metadata$assay_names <- character(0)
    }
    
    # Get the dimensions of the Seurat object
    metadata$seurat_dims <- dim(seurat_obj)
    
    # Get the number of genes
    metadata$num_genes <- metadata$seurat_dims[1]
    
    # Get the number of cells
    metadata$num_cells <- metadata$seurat_dims[2]
    
    # Get the list of dimensional reductions
    metadata$dimensional_reductions <- names(seurat_obj@reductions)
    if (is.null(metadata$dimensional_reductions)) {
        # Return an empty array (vector)
        metadata$dimensional_reductions <- character(0)
    }
    
    # Return the metadata list
    metadata
}