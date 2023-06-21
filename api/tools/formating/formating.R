library(scater)
library(anndata)
library(Seurat)
library(SingleCellExperiment)
library(SeuratDisk)
library(SeuratData)
library(patchwork)
# library(loomR)

# Get suffix of the file
get_suffix <- function(path) {
    filename <- basename(path)
    parts <- strsplit(filename,".",fixed = TRUE)
    nparts <- length(parts[[1]])
    return(parts[[1]][nparts])
}

# Read csv/xlsx/h5ad/hdf5/h5/loom/mtx/txt/tab/data/gz file to create AnnData object
load_anndata <- function(path) {
    adata <- NULL
    suffix <- tolower(get_suffix(path))
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
        delim <- detect_delim(path)
        adata <- read_text(path, delimiter = delim)
    } else if(suffix == "gz"){
        adata <- read_umi_tools(path)
    } else if(suffix == "h5Seurat" || suffix == "h5seurat"){
        Convert(path, dest = "h5ad", overwrite = TRUE, verbose = FALSE)
        adata <- read_h5ad(paste0(tools::file_path_sans_ext(path), ".h5ad"))
    } else if(suffix == "rds"){
        seurat_object <- load_seurat(path)
        SaveH5Seurat(seurat_object, overwrite = TRUE)
        Convert(paste0(tools::file_path_sans_ext(path), ".h5Seurat"), dest = "h5ad" , overwrite = TRUE, verbose = FALSE)
        adata <- read_h5ad(paste0(tools::file_path_sans_ext(path), ".h5ad"))
    } 
    # else if(suffix == "loom"){
    #     seurat_object <- load_seurat(path)
    #     SaveH5Seurat(seurat_object, overwrite = TRUE)
    #     Convert(paste0(tools::file_path_sans_ext(path), ".h5Seurat"), dest = "h5ad")
    #     ad <- read_h5ad(paste0(tools::file_path_sans_ext(path), ".h5ad"))
    # } 
    py_to_r_ifneedbe(adata)
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
        print("Inside load_seurat")
        seurat_object <- LoadH5Seurat(path)
    } else if(suffix == "h5ad"){
        Convert(path, "h5seurat", overwrite = TRUE, assay = "RNA")
        seurat_object <- LoadH5Seurat(paste0(tools::file_path_sans_ext(path), ".h5Seurat"))
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

load_sce <- function(path){
    sce <- NULL
    if(file_test("-d", path)){
        print("Inside load_seurat 3")
        if(file.exists(file.path(path,"molecules.txt")) && file.exists(file.path(path,"annotation.txt"))){
            delim <- detect_delim(path)
            molecules <- read.delim(file.path(path,"molecules.txt"), sep = delim, row.names = 1) 
            annotation <- read.delim(file.path(path,"annotation.txt"), sep = delim, stringsAsFactors = T)
            sce <- SingleCellExperiment(assays = list(counts = as.matrix(molecules)), colData = annotation)
        } else{
            print("Inside load_seurat 4")
            seurat_object <- load_seurat(path)
            if(!is.null(seurat_object)){
                sce <- as.SingleCellExperiment(seurat_object)
                rm(seurat_object) # Erase expression_matrix from memory to save RAM
            }
        }
    } else{
        print("Inside load_seurat 5")
        suffix <- tolower(get_suffix(path))
        if(suffix == "rds"){
            sce <- readRDS(path)
            print("Inside load_seurat 6")
        } else{
            print("Inside load_seurat 7")
            seurat_object <- load_seurat(path)
            if(!is.null(seurat_object)){
                sce <- as.SingleCellExperiment(seurat_object)
                rm(seurat_object) # Erase expression_matrix from memory to save RAM
            }
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
}

py_to_r_ifneedbe <- function(x) {
    if (inherits(x, "python.builtin.object")) {
        py_to_r(x)
    } else {
        x
    }
}

seurat_to_csv <- function(seurat_object, srat_path, assay = 'RNA', slot = "counts"){
    if(assay != 'RNA') slot ="data"
    csv_path <- gsub(".h5Seurat", paste("_", assay, ".csv", sep = ""), srat_path)
    if(!file.exists(csv_path)){
        write.table(as.matrix(GetAssayData(object = seurat_object, assay = assay, slot = slot)), 
        csv_path, sep = ',', row.names = T, col.names = T, quote = F)
    } else{
        print("CSV file already exists.")
    }
    csv_path
}