library(DuoClustering2018)
library(cidr)
library(Seurat)
library(SingleCellExperiment)
library(optparse)
library(dplyr)


args <- commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("-d", "--dataset"), type="character", default=NULL,
              help="input dataset name (parts before .rds)", metavar="character"),
  make_option(c("-m", "--method"), type="character", default=NULL,
              help="method name", metavar="character"),
  make_option(c("-w", "--workspace"), type="character", default="./",
              help="work directory", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list,add_help_option=FALSE)
opt = parse_args(opt_parser)

if (is.null(opt$dataset)){
  print_help(opt_parser)
  stop("Evaluation dataset must be supplied (-i).", call.=FALSE)
}
if (is.null(opt$method)){
  print_help(opt_parser)
  stop("Evaluation method must be supplied (-i).", call.=FALSE)
}



apply_Seurat <- function(sce, params, resolution) {
  (seed <- round(1e6*runif(1)))
  tryCatch({
    dat <- counts(sce)
    st <- system.time({
      data <- CreateSeuratObject(dat, min.cells = params$min.cells,
                                 min.genes = params$min.genes, project = "scRNAseq", 
                                 display.progress = FALSE) 
      data <- NormalizeData(object = data, normalization.method = "LogNormalize", 
                            scale.factor = 1e4, display.progress = FALSE)
      data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
      data <- ScaleData(object = data, display.progress = FALSE)
      data <- RunPCA(object = data, pc.genes = rownames(data@data), do.print = FALSE, 
                     pcs.compute = max(params$dims.use), seed.use = seed)
      data <- FindNeighbors(data, dims = params$dims.use)
      data <- FindClusters(data, resolution = resolution)
      cluster <- data$seurat_clusters
    })
    
    st <- c(user.self = st[["user.self"]], sys.self = st[["sys.self"]], 
            user.child = st[["user.child"]], sys.child = st[["sys.child"]],
            elapsed = st[["elapsed"]])
    list(st = st, cluster = cluster, est_k = NA)
  }, error = function(e) {
    list(st = c(user.self = NA, sys.self = NA, user.child = NA, sys.child = NA,
                elapsed = NA), 
         cluster = structure(rep(NA, ncol(sce)), names = colnames(sce)),
         est_k = NA)
  })
}


apply_CIDR <- function(sce, params, k) {
  tryCatch({
    dat <- counts(sce)
    st <- system.time({
      sData <- scDataConstructor(dat, tagType = "raw")
      sData <- determineDropoutCandidates(sData)
      sData <- wThreshold(sData)
      sData <- scDissim(sData, threads = 1)
      sData <- scPCA(sData, plotPC = FALSE)
      sData <- nPC(sData)
      
      ## Cluster with preset number of clusters
      sDataC <- scCluster(object = sData, nCluster = k, 
                          nPC = sData@nPC, cMethod = "ward.D2")
      cluster <- sDataC@clusters
      names(cluster) <- colnames(sDataC@tags)
    })
    ## Determine number of clusters automatically
    sDataA <- scCluster(object = sData, n = max(params$range_clusters),
                        nPC = sData@nPC, cMethod = "ward.D2")
    est_k <- sDataA@nCluster
    
    st <- c(user.self = st[["user.self"]], sys.self = st[["sys.self"]], 
            user.child = st[["user.child"]], sys.child = st[["sys.child"]],
            elapsed = st[["elapsed"]])
    list(st = st, cluster = cluster, est_k = est_k)
  }, error = function(e) {
    list(st = c(user.self = NA, sys.self = NA, user.child = NA, sys.child = NA,
                elapsed = NA), 
         cluster = structure(rep(NA, ncol(sce)), names = colnames(sce)),
         est_k = NA)
  })
}



method <- opt$method
scename <- opt$dataset
mydata_sce=readRDS(paste0(opt$dataset,".rds"))
params <- duo_clustering_all_parameter_settings_v2()[[paste0("sce_filteredExpr10_Koh_", 
                                                             method)]]

n_rep <- 2

set.seed(1234)
L <- lapply(seq_len(n_rep), function(i) {  ## For each run
  cat(paste0("run = ", i, "\n"))
  if (method == "Seurat") {
    tmp <- lapply(params$range_resolutions, function(resolution) {  
      ## For each resolution
      cat(paste0("resolution = ", resolution, "\n"))
      ## Run clustering
      res <- get(paste0("apply_", method))(sce = mydata_sce, params = params, 
                                           resolution = resolution)
      
      ## Put output in data frame
      df <- data.frame(dataset = scename, 
                       method = method, 
                       cell = names(res$cluster),
                       run = i,
                       k = length(unique(res$cluster)),
                       resolution = resolution,
                       cluster = res$cluster,
                       stringsAsFactors = FALSE, row.names = NULL)
      tm <- data.frame(dataset = scename, 
                       method = method,
                       run = i, 
                       k = length(unique(res$cluster)),
                       resolution = resolution,
                       user.self = res$st[["user.self"]],
                       sys.self = res$st[["sys.self"]],
                       user.child = res$st[["user.child"]],
                       sys.child = res$st[["sys.child"]],
                       elapsed = res$st[["elapsed"]],
                       stringsAsFactors = FALSE, row.names = NULL)
      kest <- data.frame(dataset = scename, 
                         method = method,
                         run = i, 
                         k = length(unique(res$cluster)),
                         resolution = resolution,
                         est_k = res$est_k,
                         stringsAsFactors = FALSE, row.names = NULL)
      list(clusters = df, timing = tm, kest = kest)
    })  ## End for each resolution
  } else {
    tmp <- lapply(params$range_clusters, function(k) {  ## For each k
      cat(paste0("k = ", k, "\n"))
      ## Run clustering
      res <- get(paste0("apply_", method))(sce = mydata_sce, params = params, k = k)
      
      ## Put output in data frame
      df <- data.frame(dataset = scename, 
                       method = method, 
                       cell = names(res$cluster),
                       run = i,
                       k = k,
                       resolution = NA,
                       cluster = res$cluster,
                       stringsAsFactors = FALSE, row.names = NULL)
      tm <- data.frame(dataset = scename, 
                       method = method,
                       run = i, 
                       k = k,
                       resolution = NA,
                       user.self = res$st[["user.self"]],
                       sys.self = res$st[["sys.self"]],
                       user.child = res$st[["user.child"]],
                       sys.child = res$st[["sys.child"]],
                       elapsed = res$st[["elapsed"]],
                       stringsAsFactors = FALSE, row.names = NULL)
      kest <- data.frame(dataset = scename, 
                         method = method,
                         run = i, 
                         k = k,
                         resolution = NA,
                         est_k = res$est_k,
                         stringsAsFactors = FALSE, row.names = NULL)
      list(clusters = df, timing = tm, kest = kest)
    })  ## End for each k
  }
  
  ## Summarize across different values of k
  assignments <- do.call(rbind, lapply(tmp, function(w) w$clusters))
  timings <- do.call(rbind, lapply(tmp, function(w) w$timing))
  k_estimates <- do.call(rbind, lapply(tmp, function(w) w$kest))
  list(assignments = assignments, timings = timings, k_estimates = k_estimates)
})  ## End for each run

assignments <- do.call(rbind, lapply(L, function(w) w$assignments))
timings <- do.call(rbind, lapply(L, function(w) w$timings))
k_estimates <- do.call(rbind, lapply(L, function(w) w$k_estimates))

assignments$trueclass <- mydata_sce$Tru_lab$majorType[match(assignments$cell, mydata_sce$Tru_lab$X)]

res <- list(assignments = assignments, timings = timings,
            k_estimates = k_estimates)

df <- dplyr::full_join(res$assignments %>%
                         dplyr::select(dataset, method, cell, run, k, 
                                       resolution, cluster, trueclass),
                       res$k_estimates %>%
                         dplyr::select(dataset, method, run, k, 
                                       resolution, est_k)
) %>% dplyr::full_join(res$timings %>% dplyr::select(dataset, method, run, k,
                                                     resolution, elapsed))

setwd(opt$workspace)

result=df
rm(df)
if(file.exists("df.RData")){
  load("df.RData")
  df=rbind(df,result)
}else{
  df=result
}
save(df,file="df.RData")
print("Done")