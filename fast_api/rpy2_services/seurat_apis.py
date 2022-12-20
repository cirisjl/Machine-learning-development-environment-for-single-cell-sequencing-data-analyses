import os
from fastapi import APIRouter
from rpy2 import robjects

router = APIRouter(
    prefix="/rpy2_services",
    tags=["rpy2_services"],
    # dependencies=[Depends(jwt.JWTBearer())],
    responses={404: {"description": "Error in calling rpy2 service"}},
)

@router.get("/qc_seurat")
def qc_seurat(input_path, dataset_name, output_path, source):
    if os.path.exists(input_path):
        input_filename = input_path+dataset_name
        if os.path.exists(input_filename):
            robjects.r(
				'''
                        qc_f <- function(input_path,input_filename, output_path,src) {
                            library(Seurat)
							library(ggplot2)
							library(SingleR)
							library(dplyr)
							library(celldex)
							library(RColorBrewer)
							library(SingleCellExperiment)
							library(SeuratDisk)
							library(SeuratData)
							library(anndata)
							library(patchwork)
							source(src)

							input= input_path
							output_seurat= output_path
							srat <- load_seurat(input)
							str(srat)

							meta <- srat@meta.data
							dim(meta)

							head(meta)

							summary(meta$nCount_RNA)

							summary(meta$nFeature_RNA)

							srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^MT-")

							srat[["percent.rb"]] <- PercentageFeatureSet(srat, pattern = "^RP[SL]")


							doublets <- read.table(input_filename,header = F,row.names = 1)
							colnames(doublets) <- c("Doublet_score","Is_doublet")
							srat <- AddMetaData(srat,doublets)
							head(srat[[]])

							png('rplot.png')
							VlnPlot(srat, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"),ncol = 4,pt.size = 0.1) &
							theme(plot.title = element_text(size=10))
							dev.off()

							FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "percent.mt")

							FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

							FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "percent.rb")

							FeatureScatter(srat, feature1 = "percent.rb", feature2 = "percent.mt")

							FeatureScatter(srat, feature1 = "nFeature_RNA",
							               feature2 = "Doublet_score")

							srat[['QC']] <- ifelse(srat@meta.data$Is_doublet == 'True','Doublet','Pass')
							srat[['QC']] <- ifelse(srat@meta.data$nFeature_RNA < 500 & srat@meta.data$QC == 'Pass','Low_nFeature',srat@meta.data$QC)
							srat[['QC']] <- ifelse(srat@meta.data$nFeature_RNA < 500 & srat@meta.data$QC != 'Pass' & srat@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',srat@meta.data$QC,sep = ','),srat@meta.data$QC)
							srat[['QC']] <- ifelse(srat@meta.data$percent.mt > 15 & srat@meta.data$QC == 'Pass','High_MT',srat@meta.data$QC)
							srat[['QC']] <- ifelse(srat@meta.data$nFeature_RNA < 500 & srat@meta.data$QC != 'Pass' & srat@meta.data$QC != 'High_MT',paste('High_MT',srat@meta.data$QC,sep = ','),srat@meta.data$QC)
							table(srat[['QC']])

							VlnPlot(subset(srat, subset = QC == 'Pass'),
									features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 4, pt.size = 0.1) &
							theme(plot.title = element_text(size=10))


							#SaveH5Seurat(srat, filename = output_seurat, overwrite = TRUE)
							# Convert(srat, dest = "h5ad", overwrite = TRUE)


                        }
                        '''
						)
            r_qc_f = robjects.globalenv['qc_f']
            qc_res = r_qc_f(input_path, input_filename, output_path, source)
            print(qc_res)
			
        else:
			#print("Dataset given doesn't exists")
            raise NameError("Dataset given doesn't exists")  # Raise Error
    else:
        #print(" Input Directory doesn't exists")
        raise NameError("Input Directory doesn't exists")  # Raise Error


if __name__ == '__main__':
    ip = "D:/Master's/ML_backend_AISingleCell/Machine-learning-development-environment-for-single-cell-sequencing-data-analyses/fast_api/Seurat_data/input_data"
    op = "D:/Master's/ML_backend_AISingleCell/Machine-learning-development-environment-for-single-cell-sequencing-data-analyses/fast_api/Seurat_data/output_seurat/output_seurat"
    dn = "/scrublet_calls.tsv"
    src = "D:/Master's/ML_backend_AISingleCell/Machine-learning-development-environment-for-single-cell-sequencing-data-analyses/fast_api/Seurat_data/formatting.R"
    res = qc_seurat(ip, dn, op, src)
    print(res)
