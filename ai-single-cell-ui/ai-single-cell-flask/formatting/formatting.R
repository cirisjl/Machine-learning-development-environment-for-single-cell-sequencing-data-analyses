library(scater)
library(anndata)
library(Seurat)
library(SingleCellExperiment)
library(SeuratDisk)
library(SeuratData)
library(patchwork)
library(Signac)

# Get suffix of the file
get_suffix <- function(path) {
    filename <- basename(path)
    parts <- strsplit(filename,".",fixed = TRUE)
    nparts <- length(parts[[1]])
    return(parts[[1]][nparts])
}

load_expression_matrix <- function(path){
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
            delim <- detect_delim(path)
            molecules <- read.delim(file.path(path,"molecules.txt"), sep = delim, row.names = 1) 
            expression_matrix <- as.matrix(molecules)} 
    } else{       
        suffix <- tolower(get_suffix(path))
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

load_seurat <- function(path, project = NULL){
    seurat_object <- NULL
    suffix <- tolower(get_suffix(path))
    print("Inside load_seurat 1")

    if(suffix == "h5Seurat" || suffix == "h5seurat"){
        seurat_object <- LoadH5Seurat(path)
    } else if(suffix == "h5ad"){
        Convert(path, "h5seurat", overwrite = TRUE, assay = "RNA")
        seurat_object <- LoadH5Seurat(paste0(tools::file_path_sans_ext(path), ".h5seurat"))
    } else if(suffix == "rds"){
        sce <- readRDS(path)
        seurat_object <- as.Seurat(sce, slot = "counts", data = NULL)
    } else {
        expression_matrix <- load_expression_matrix(path)
        if(!is.null(expression_matrix) && !is.null(project)) {
            seurat_object <- CreateSeuratObject(counts = expression_matrix, project = project)
        } else if (!is.null(expression_matrix)){
            seurat_object <- CreateSeuratObject(counts = expression_matrix)
        }
        rm(expression_matrix) # Erase expression_matrix from memory to save RAM
    }
    # else if(suffix == "loom"){
    #     loom <- connect(filename = path, mode = "r")
    #     seurat_object <- as.Seurat(loom)
    # }
    seurat_object
}

load_metadata <- function(seurat_obj) {
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
detect_delim <- function(path, nchar = 1e3) {
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

convert_to_anndata <- function(path, assay = 'RNA') {
    adata_path <- NULL
    suffix <- tolower(get_suffix(path))
    if(suffix == "h5Seurat" || suffix == "h5seurat"){
        if(assay != 'RNA') {
            DefaultAssay(seurat_object) <- assay
            SaveH5Seurat(seurat_object, filename = path, overwrite = TRUE, verbose = FALSE)
        }
        adata_path <- Convert(path, dest = "h5ad", assay=assay, overwrite = TRUE, verbose = FALSE)
    } else if(suffix == "rds"){
        seurat_object <- load_seurat(path)
        seurat_path <- paste0(tools::file_path_sans_ext(path), ".h5Seurat")
        SaveH5Seurat(seurat_object, filename = seurat_path, overwrite = TRUE, verbose = FALSE)
        adata_path <- Convert(seurat_path, dest = "h5ad" , overwrite = TRUE, verbose = FALSE)
    } 
    adata_path
}



# Return a list of assays if the default assay is not "RNA"
convert_seurat_sce_to_anndata <- function(path, assay = NULL) {
    seurat_object <- NULL
    assay_names <- NULL
    default_assay <- NULL
    anndata_path <- NULL
    suffix <- tolower(get_suffix(path))

    if(suffix == "h5Seurat" || suffix == "h5seurat"){ # Seurat object
        seurat_object <- LoadH5Seurat(path)
        default_assay <- DefaultAssay(seurat_object)
        assay_names <- names(seurat_object@assays)

        # Return a list of assays if the default assay is not "RNA"
        if (default_assay == 'RNA') {
            anndata_path <- convert_to_anndata(path)
            # anndata_path <- save_as_anndata(seurat_object, path, assay = 'RNA')
        } else if (!is.null(assay) && (assay %in% assay_names)){
            # anndata_path <- save_as_anndata(seurat_object, path, assay = assay)
            anndata_path <- convert_to_anndata(path, assay = assay)
        } 
        seurat_object <- NULL
    } else if (suffix == "rds") { # Single-Cell Experiment object
        anndata_path <- convert_to_anndata(path)
    }

    list(default_assay=default_assay, assay_names=assay_names, anndata_path=anndata_path) # First, check if the anndata_path is NULL; Second, check the assay_names
}