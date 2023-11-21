# Read matrix file to expFile
# expFile is the count matrix, where rownames are gene ID or symbols, column names are cells
expFile

# Options
species_name <- "Human"
species_name <- "Mouse"

suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(org.Mm.eg.db))

# detect rownames gene list with identifer by the largest number of matches: 1) gene symbol. 2)ensembl geneid. 3) ncbi entrez id. 4)
get_rowname_type <- function (l, db) {
  res1 <-
    tryCatch(
      nrow(
        AnnotationDbi::select(
          db,
          keys = l,
          columns = c("ENTREZID", "SYMBOL", "ENSEMBL", "ENSEMBLTRANS"),
          keytype = "SYMBOL"
        )
      ),
      error = function(e)
        0
    )
  res2 <-
    tryCatch(
      nrow(
        AnnotationDbi::select(
          db,
          keys = l,
          columns = c("ENTREZID", "SYMBOL", "ENSEMBL", "ENSEMBLTRANS"),
          keytype = "ENSEMBL"
        )
      ),
      error = function(e)
        0
    )
  res3 <-
    tryCatch(
      nrow(
        AnnotationDbi::select(
          db,
          keys = l,
          columns = c("ENTREZID", "SYMBOL", "ENSEMBL", "ENSEMBLTRANS"),
          keytype = "ENTREZID"
        )
      ),
      error = function(e)
        0
    )
  res4 <-
    tryCatch(
      nrow(
        AnnotationDbi::select(
          db,
          keys = l,
          columns = c("ENTREZID", "SYMBOL", "ENSEMBL", "ENSEMBLTRANS"),
          keytype = "ENSEMBLTRANS"
        )
      ),
      error = function(e)
        0
    )
  result <- c("error", "SYMBOL", "ENSEMBL", "ENTREZID", "ENSEMBLTRANS")
  result_vec <- c(1, res1, res2, res3, res4)
  return(c(result[which.max(result_vec)], result_vec[which.max(result_vec)]))
  #write("No matched gene identifier found, please check your data.",file=paste(jobid,"_error.txt",sep=""),append=TRUE)
}


db <- c("Mouse" = org.Mm.eg.db, "Human" = org.Hs.eg.db)
select_db <- db[which(names(db) %in% species_name)]
gene_identifier <-
  sapply(select_db, get_rowname_type, l = rownames(expFile))
main_species <- names(which.max(gene_identifier[2,]))
main_db <- db[which(names(db) %in% main_species)][[1]]
main_identifier <-
  as.character(gene_identifier[1, which.max(gene_identifier[2,])])

all_match <-
  AnnotationDbi::select(
    main_db,
    keys = rownames(expFile),
    columns = c("SYMBOL", "ENSEMBL"),
    keytype = main_identifier
  )

expFile <- merge(expFile, all_match, by.x = 0, by.y = main_identifier)
dim(expFile)
expFile <- na.omit(expFile)

## merge expression values with same gene names
if (main_identifier == "ENSEMBL") {
  expFile <- expFile[, -1]
  expFile <- aggregate(. ~ SYMBOL, expFile, sum)
} else if (main_identifier == "ENSEMBLTRANS") {
  expFile <- expFile[, c(-1, -(ncol(expFile)))]
  expFile <- aggregate(. ~ SYMBOL, expFile, sum)
} else {
  expFile <- expFile[, -(ncol(expFile))]
  expFile <- aggregate(. ~ Row.names, expFile, sum)
}
expFile <- expFile[!duplicated(expFile[, 1]), ]
rownames(expFile) <- expFile[, 1]
expFile <- expFile[, -1]

## remove rows with empty gene name
if (length(which(rownames(expFile) == "")) > 0) {
  expFile <- expFile[-which(rownames(expFile) == ""), ]
}
