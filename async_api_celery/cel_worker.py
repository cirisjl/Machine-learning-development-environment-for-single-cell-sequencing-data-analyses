from celery import Celery
import rpy2
import rpy2.robjects
import os
from rpy2.robjects.packages import importr

scater = importr('scater')
singleCellExperiment = importr('SingleCellExperiment')
seurat = importr('Seurat')
seuratDisk = importr('SeuratDisk')
patchwork = importr('patchwork')
anndata = importr('anndata')
EnsDb_Hsapiens_v86 = importr('EnsDb.Hsapiens.v86')
seuratData = importr('SeuratData')
scales = importr('scales')


os.environ.setdefault("FORKED_BY_MULTIPROCESSING",'1')

app = Celery('Celery worker', broker="amqp://localhost")

@app.task(name="start")
def start(path):
    # Use a breakpoint in the code line below to debug your script.
    rpy2.robjects.r('''
        load_sce <- function(path){
            sce <- NULL 
            if(file_test("-d", path)){
                if(file.exists(file.path(path,"molecules.txt")) && file.exists(file.path(path,"annotations.txt"))){
                    deli <- NULL
                    nchar <- 1e3
                    if (file.exists(path)) {
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
                          deli <- substr(chars, search, search + attr(search, "match.length") - 1)
                        }
                      }
                      "\r\n"
                      molecules <- read.delim(file.path(path,"molecules.txt"), sep = deli, row.names = 1)
                      annotation <- read.delim(file.path(path,"annotations.txt"), sep = deli, stringsAsFactors = T)
                      sce <- SingleCellExperiment(assays = list(counts = as.matrix(molecules)), colData = annotation)
                    }
                }else{
                    suffix <- NULL
                    filename <- basename(path)
                    parts<-strsplit(filename,".",fixed = TRUE)
                    nparts<-length(parts[[1]])
                    suffix <- parts[[1]][nparts]
                    if(suffix == "rds"){
                        expression_matrix <- readRDS(path)
                        sce <- Convert(from = expression_matrix, to = "sce")
                    }
                    print("end of else block")
                }
                print(sce)
                altExp(sce,"ERCC") <- sce[grep("^ERCC-",rownames(sce)), ]
                sce <- sce[grep("^ERCC-",rownames(sce),invert = T), ]
                ensdb_genes <- genes(EnsDb.Hsapiens.v86)
                MT_names <- ensdb_genes[seqnames(ensdb_genes) == "MT"]$gene_id
                is_mito <- rownames(sce) %in% MT_names
                table(is_mito)
                sce_cell <- perCellQCMetrics(sce,subsets=list(Mito=is_mito))
                sce_feature <- perFeatureQCMetrics(sce)
                print("Printing sce cells")
                print(head(sce_cell))
                print("Printing sce feature")
                print(head(sce_feature))
                sce <- addPerCellQC(sce, subsets=list(Mito=is_mito))
                sce <- addPerFeatureQC(sce)
                hist(
                    sce$total,
                    breaks = 100
                )
                abline(v = 25000, col = "red")

                hist(
                  sce_cell$detected,
                  breaks = 100
                )
                abline(v = 7000, col = "red")
                reasons <- quickPerCellQC(sce_cell, sub.fields=c("subsets_Mito_percent", "altexps_ERCC_percent"))
                colSums(as.matrix(reasons))
                sce$discard <- reasons$discard
                plotColData(sce, x="sum", y="subsets_Mito_percent", colour_by="discard")
                plotColData(sce, x="sum", y="detected", colour_by="discard")
                plotColData(sce, x="altexps_ERCC_percent", y="subsets_Mito_percent",colour_by="discard")
                plotColData(sce, x="sum", y="detected", colour_by="discard", other_fields = "individual") + 
                facet_wrap(~individual) + scale_x_continuous(labels = unit_format(unit = "k", scale = 1e-3))
                plotColData(sce, x="sum", y="detected", colour_by="discard", other_fields = "replicate") + 
                facet_wrap(~replicate)  + scale_x_continuous(labels = unit_format(unit = "k", scale = 1e-3))

                plotHighestExprs(sce, exprs_values = "counts", feature_names_to_plot = NULL, colour_cells_by="detected")
                keep_feature <- nexprs(sce, byrow = TRUE, detection_limit = 1) >= 2
                rowData(sce)$discard <- ! keep_feature
                table(rowData(sce)$discard)
                sce.qc <- sce[! rowData(sce)$discard,! colData(sce)$discard]
                srat <- as.Seurat(sce.qc, counts = "counts", data = "counts")
                srat <- RenameAssays(object = srat, originalexp = 'RNA')
                DefaultAssay(object = srat) <- "RNA"
                SaveH5Seurat(srat, filename = "output_seurat", overwrite = TRUE, file=file)

                srat

        }

    ''')
    function_call = rpy2.robjects.r['load_sce']
    sert = function_call(path)

    print("Seurat Object : ")
    print(sert)
    print("completed")
    return {"message": "success"}

@app.task(name="test")
def test(temp):
    print("temp is ")
    return temp[::-1]

