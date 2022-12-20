library(scater)
library(anndata)
library(Seurat)
library(SingleCellExperiment)
library(SeuratDisk)
library(SeuratData)
library(patchwork)
library(loomR)

# Get suffix of the file
get_suffix <- function(path) {
    filename <- basename(path)
    parts<-strsplit(filename,".",fixed = TRUE)
    nparts<-length(parts[[1]])
    return(parts[[1]][nparts])
}

# Read csv/xlsx/h5ad/hdf5/h5/loom/mtx/txt/tab/data/gz file to create AnnData object
load_anndata <- function(path) {
    ad <- NULL
    suffix <- get_suffix(path)
    if(suffix == "csv"){
        ad <- read_csv(path)
    } else if(suffix == "xlsx"){
        ad <- read_excel(path)
    } else if(suffix == "h5ad"){
        ad <- read_h5ad(path)
    } else if(suffix == "hdf5" || suffix == "h5"){
        ad <- read_hdf(path)
    } else if(suffix == "loom"){
        ad <- read_loom(path)
    } else if(suffix == "mtx"){
        ad <- read_mtx(path)
    } else if(suffix == "txt" || suffix == "tab" || suffix == "data"){
        delim <- detect_delim(path)
        ad <- read_text(path, delimiter = delim)
    } else if(suffix == "gz"){
        ad <- read_umi_tools(path)
    } else if(suffix == "h5Seurat" || suffix == "h5seurat"){
        Convert(path, dest = "h5ad")
        ad <- read_h5ad(paste0(tools::file_path_sans_ext(path), ".h5ad"))
    } else if(suffix == "rds"){
        seurat_object <- load_seurat(path)
        SaveH5Seurat(seurat_object, overwrite = TRUE)
        Convert(paste0(tools::file_path_sans_ext(path), ".h5Seurat"), dest = "h5ad")
        ad <- read_h5ad(paste0(tools::file_path_sans_ext(path), ".h5ad"))
    } else if(suffix == "loom"){
        seurat_object <- load_seurat(path)
        SaveH5Seurat(seurat_object, overwrite = TRUE)
        Convert(paste0(tools::file_path_sans_ext(path), ".h5Seurat"), dest = "h5ad")
        ad <- read_h5ad(paste0(tools::file_path_sans_ext(path), ".h5ad"))
    } 
    ad
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
        suffix <- get_suffix(path)
        if(suffix == "h5"){
            expression_matrix  <- Read10X_h5(path)
        }
        else if(suffix == "csv"){
            expression_matrix  <- as.matrix(read.csv(path, header=TRUE, row.names=1))
        } 
    }
    expression_matrix
}

load_seurat <- function(path){
    seurat_object <- NULL
    expression_matrix <- load_expression_matrix(path)
    if(is.null(expression_matrix)){
        suffix <- get_suffix(path)
        if(suffix == "h5Seurat" || suffix == "h5seurat"){
            seurat_object <- LoadH5Seurat(path)
        } else if(suffix == "h5ad"){
            Convert(path, "h5seurat", overwrite = TRUE, assay = "RNA")
            seurat_object <- LoadH5Seurat(paste0(tools::file_path_sans_ext(path), ".h5Seurat"))
        } else if(suffix == "rds"){
            expression_matrix <- readRDS(path)
            seurat_object <- as.Seurat(expression_matrix, counts = "counts", data = "logcounts")
        } else if(suffix == "loom"){
            loom <- connect(filename = path, mode = "r")
            seurat_object <- as.Seurat(loom)
        }
    } else{
        seurat_object <- CreateSeuratObject(counts = expression_matrix)
    }
    expression_matrix <- NULL # Erase expression_matrix from memory to save RAM,
    seurat_object
}

load_sce <- function(path){
    sce <- NULL
    if(file_test("-d", path)){
        if(file.exists(file.path(path,"molecules.txt")) && file.exists(file.path(path,"annotation.txt"))){
            delim <- detect_delim(path)
            molecules <- read.delim(file.path(path,"molecules.txt"), sep = delim, row.names = 1) 
            annotation <- read.delim(ile.path(path,"annotation.txt"), sep = delim, stringsAsFactors = T)
            sce <- SingleCellExperiment(assays = list(counts = as.matrix(molecules)), colData = annotation)
        }
    } else{
        suffix <- get_suffix(path)
        if(suffix == "rds"){
            expression_matrix <- readRDS(path)
            sce <- Convert(from = expression_matrix, to = "sce")
        }
    }
    sce

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
  "\r\n"
}