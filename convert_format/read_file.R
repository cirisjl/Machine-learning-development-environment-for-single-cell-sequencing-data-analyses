#' Read file switcher
#'
#' @param type string
#' @param path string
#'
#' @return dataframe/matrix/sparse matrix
#'

read_data <-
  function(type = "text",
           # originalname = 'pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5',
           path = "c57e41d078ce9810717de392b4451605.h5ad") {
    path_prefix <- "/data/"
    absolute_path <- paste0(path_prefix, path)
    result <-
      switch(type,
        "text" = read_data_text(absolute_path),
        "10x_h5" = read_data_10xh5(absolute_path),
        "h5ad" = read_data_h5ad(absolute_path)
      )
    return(result)
  }

#' Read csv,tsv,txt format
#'
#' @param path string
#'
#' @return dataframe/matrix/sparse matrix
#'

read_data_text <- function(path) {
  delim <- detect_delim(path)
  result <- read.table(path,
    sep = delim,
    header = T,
    row.names = 1
  )
  return(result)
}

#' Read 10x hdf5 format
#'
#' @param path string
#'
#' @return dataframe/matrix/sparse matrix
#'

read_data_10xh5 <- function(path) {
  result <- Read10X_h5(path)
  return(result)
}

#' Read anndata object from h5ad
#'
#' @param path string
#'
#' @return dataframe/matrix/sparse matrix
#'

read_data_h5ad <- function(path) {
  require(SeuratDisk)
  require(tools)
  Convert(path, ".h5seurat")
  result <- LoadH5Seurat(paste0(tools::file_path_sans_ext(path), ".h5Seurat"))
  return(result)
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
