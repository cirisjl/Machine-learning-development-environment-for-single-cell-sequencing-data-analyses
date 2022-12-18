from celery import Celery
import rpy2
import rpy2.robjects
import os
from rpy2.robjects.packages import importr

patchwork = importr('patchwork')
harmony = importr('harmony')
seurat = importr('Seurat')
seuratDisk = importr('SeuratDisk')
SeuratWrappers = importr('SeuratWrappers')
rliger = importr('rliger')
reshape = importr('reshape2')
rColorBrewer = importr('RColorBrewer')
dplyr = importr('dplyr')

os.environ.setdefault("FORKED_BY_MULTIPROCESSING",'1')

app = Celery('Celery worker', broker="amqp://localhost",backend='')

@app.task(name="start")
def start(dataset1, dataset2):
    rpy2.robjects.r('''
        source("F:/Lab(Single_cell_analysis)/custom_seurat_functions.R")
        source("F:/Lab(Single_cell_analysis)/formatting.R")
        batch_correction <- function(dataset1, dataset2){
            print("Dataset 1 is :")
            print(dataset1)
            print("Dataset 2 is :")
            print(dataset2)
            matrix_3p <- load_expression_matrix(dataset1)
            matrix_5p <- load_expression_matrix(dataset2)$`Gene Expression`
            srat_3p   <- CreateSeuratObject(matrix_3p,project = "pbmc10k_3p")
            srat_5p   <- CreateSeuratObject(matrix_5p,project = "pbmc10k_5p")
            rm(matrix_3p)
            rm(matrix_5p)

            srat_3p[["percent.mt"]]  <- PercentageFeatureSet(srat_3p, pattern = "^MT-")
            srat_3p[["percent.rbp"]] <- PercentageFeatureSet(srat_3p, pattern = "^RP[SL]")
            srat_5p[["percent.mt"]]  <- PercentageFeatureSet(srat_5p, pattern = "^MT-")
            srat_5p[["percent.rbp"]] <- PercentageFeatureSet(srat_5p, pattern = "^RP[SL]")

            VlnPlot(srat_3p, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)

            VlnPlot(srat_5p, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)

            table(rownames(srat_3p) %in% rownames(srat_5p)) 

            srat_3p <- subset(srat_3p, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 15)
            srat_5p <- subset(srat_5p, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)

            pbmc_list <- list()
            pbmc_list[["pbmc10k_3p"]] <- srat_3p
            pbmc_list[["pbmc10k_5p"]] <- srat_5p

            for (i in 1:length(pbmc_list)) {
              pbmc_list[[i]] <- NormalizeData(pbmc_list[[i]], verbose = F)
              pbmc_list[[i]] <- FindVariableFeatures(pbmc_list[[i]], selection.method = "vst", nfeatures = 2000, verbose = F)
            }

            pbmc_anchors    <- FindIntegrationAnchors(object.list = pbmc_list, dims = 1:30)

            pbmc_seurat     <- IntegrateData(anchorset = pbmc_anchors, dims = 1:30)

            DefaultAssay(pbmc_seurat) <- "RNA"

            pbmc_seurat <- NormalizeData(pbmc_seurat, verbose = F)
            pbmc_seurat <- FindVariableFeatures(pbmc_seurat, selection.method = "vst", nfeatures = 2000, verbose = F)
            pbmc_seurat <- ScaleData(pbmc_seurat, verbose = F)
            pbmc_seurat <- RunPCA(pbmc_seurat, npcs = 30, verbose = F)
            pbmc_seurat <- RunUMAP(pbmc_seurat, reduction = "pca", dims = 1:30, verbose = F)

            DimPlot(pbmc_seurat,reduction = "umap") + plot_annotation(title = "Seurat 3 integration",
                                                                      theme = theme(plot.title = element_text(hjust = 0.5)))

            DefaultAssay(pbmc_seurat) <- "integrated"
            pbmc_seurat <- ScaleData(pbmc_seurat, verbose = F)
            pbmc_seurat <- RunPCA(pbmc_seurat, npcs = 30, verbose = F)
            pbmc_seurat <- RunUMAP(pbmc_seurat, reduction = "pca", dims = 1:30, verbose = F)

            DimPlot(pbmc_seurat, reduction = "umap") + plot_annotation(title = "Seurat 3 integration",
                                                           theme = theme(plot.title = element_text(hjust = 0.5)))

            DimPlot(pbmc_seurat, reduction = "umap", split.by = "orig.ident") + NoLegend()

            pbmc_seurat <- FindNeighbors(pbmc_seurat, dims = 1:30, k.param = 10, verbose = F)
            pbmc_seurat <- FindClusters(pbmc_seurat, verbose = F)
            DimPlot(pbmc_seurat,label = T) + NoLegend()

            count_table <- table(pbmc_seurat@meta.data$seurat_clusters, pbmc_seurat@meta.data$orig.ident)
            count_table

            plot_integrated_clusters(pbmc_seurat) 

            #Suerat V4

            features <- SelectIntegrationFeatures(object.list = pbmc_list)

            immune.anchors <- FindIntegrationAnchors(object.list = pbmc_list, anchor.features = features)

            immune.combined <- IntegrateData(anchorset = immune.anchors)

            DefaultAssay(immune.combined) <- "integrated"

            # Run the standard workflow for visualization and clustering

            immune.combined <- ScaleData(immune.combined, verbose = FALSE)
            immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
            immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
            immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
            immune.combined <- FindClusters(immune.combined, resolution = 0.5)

            immune.combined    <- SetIdent(immune.combined,value = "orig.ident")

            DimPlot(immune.combined, reduction = "umap")+ plot_annotation(title = "Seurat 4 integration",
                                                                          theme = theme(plot.title = element_text(hjust = 0.5)))


            DimPlot(immune.combined, reduction = "umap", label = FALSE, repel = TRUE)+ plot_annotation(title = "Seurat integration",
                                                                                           theme = theme(plot.title = element_text(hjust = 0.5)))

            #Harmony, 3’ vs 5’ 10k PBMC

            pbmc_harmony    <- merge(srat_3p,srat_5p)

            pbmc_harmony <- NormalizeData(pbmc_harmony, verbose = F)
            pbmc_harmony <- FindVariableFeatures(pbmc_harmony, selection.method = "vst", nfeatures = 2000, verbose = F)
            pbmc_harmony <- ScaleData(pbmc_harmony, verbose = F)
            pbmc_harmony <- RunPCA(pbmc_harmony, npcs = 30, verbose = F)
            pbmc_harmony <- RunUMAP(pbmc_harmony, reduction = "pca", dims = 1:30, verbose = F)

            DimPlot(pbmc_harmony,reduction = "umap") + plot_annotation(title = "Before integration",
                                                           theme = theme(plot.title = element_text(hjust = 0.5)))

            pbmc_harmony <- pbmc_harmony %>% RunHarmony("orig.ident", plot_convergence = T)

            harmony_embeddings <- Embeddings(pbmc_harmony, 'harmony')
            harmony_embeddings[1:5, 1:5]

            p1 <- DimPlot(object = pbmc_harmony, reduction = "harmony", pt.size = .1, group.by = "orig.ident") + NoLegend()
            p2 <- VlnPlot(object = pbmc_harmony, features = "harmony_1", group.by = "orig.ident", pt.size = .1) + NoLegend()
            plot_grid(p1,p2)

            pbmc_harmony <- pbmc_harmony %>% 
              RunUMAP(reduction = "harmony", dims = 1:30, verbose = F) %>% 
              FindNeighbors(reduction = "harmony", k.param = 10, dims = 1:30) %>% 
              FindClusters() %>% 
              identity()

            pbmc_harmony <- SetIdent(pbmc_harmony,value = "orig.ident")
            DimPlot(pbmc_harmony,reduction = "umap") + plot_annotation(title = "Harmony integration",
                                                                       theme = theme(plot.title = element_text(hjust = 0.5)))

            DimPlot(pbmc_harmony, reduction = "umap", group.by = "orig.ident", pt.size = .1, split.by = 'orig.ident') + NoLegend()

            pbmc_harmony <- SetIdent(pbmc_harmony,value = "seurat_clusters")
            DimPlot(pbmc_harmony,label = T) + NoLegend()

            plot_integrated_clusters(pbmc_harmony)

            rm(pbmc_harmony)

            #LIGER, 3’ vs 5’ 10k PBMC

            pbmc_liger    <- merge(srat_3p,srat_5p)

            pbmc_liger    <- NormalizeData(pbmc_liger)
            pbmc_liger    <- FindVariableFeatures(pbmc_liger)
            pbmc_liger    <- ScaleData(pbmc_liger, split.by = "orig.ident", do.center = F)

            pbmc_liger    <- RunOptimizeALS(pbmc_liger, k = 30, lambda = 5, split.by = "orig.ident") ## this one takes a while

            pbmc_liger    <- RunQuantileNorm(pbmc_liger, split.by = "orig.ident")

            pbmc_liger    <- FindNeighbors(pbmc_liger,reduction = "iNMF",k.param = 10,dims = 1:30)

            pbmc_liger    <- FindClusters(pbmc_liger)

            pbmc_liger    <- RunUMAP(pbmc_liger, dims = 1:ncol(pbmc_liger[["iNMF"]]), reduction = "iNMF", verbose = F)
            pbmc_liger    <- SetIdent(pbmc_liger,value = "orig.ident")
            DimPlot(pbmc_liger,reduction = "umap") + plot_annotation(title = "LIGER integration",
                                                                     theme = theme(plot.title = element_text(hjust = 0.5)))

            DimPlot(pbmc_liger, reduction = "umap", group.by = "orig.ident", pt.size = .1, split.by = 'orig.ident') + NoLegend()

            pbmc_liger <- SetIdent(pbmc_liger,value = "seurat_clusters")
            DimPlot(pbmc_liger,reduction = "umap",label = T) + NoLegend()

            plot_integrated_clusters(pbmc_liger)

            rm(pbmc_liger)
            rm(srat_3p)
            rm(srat_5p)

        }


    ''')
    function_call = rpy2.robjects.r['batch_correction']

    function_call(dataset1, dataset2)

    print("Process Completed")
