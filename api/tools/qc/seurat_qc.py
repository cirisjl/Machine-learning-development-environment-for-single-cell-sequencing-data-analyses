import os
import scanpy as sc
import numpy as np
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from tools.formating.formating import change_file_extension, load_anndata, convert_from_r, convert_to_r
from utils.redislogger import redislogger

# Ensure that pandas2ri is activated for automatic conversion
pandas2ri.activate()

def run_seurat_qc(input, unique_id, output, assay='RNA', min_genes=200, max_genes=0, min_UMI_count=0, max_UMI_count=0, percent_mt_max=5, percent_rb_min=0, resolution=0.5, dims=10, n_pcs=3, doublet_rate=0.075, regress_cell_cycle=False):
    redislogger.info(unique_id, "Start Seurat Quality Control...")

    adata = None
    default_assay = None
    assay_names = None
    ddl_assay_names = False

    if assay is None:
        redislogger.info(unique_id, "Assay is not provides, setting it to 'RNA'.")
        assay = 'RNA'
    if max_genes is None:
        max_genes = 0

    adata_path = change_file_extension(output, 'h5ad')

    try:
        # Defining the R script and loading the instance in Python
        ro.r['source'](os.path.abspath(os.path.join(os.path.dirname(__file__), 'seurat_qc.R')))
        RunSeuratQC_r = ro.globalenv['RunSeuratQC']
        redislogger.info(unique_id, "Running R script for Seurat QC.")
        results = list(RunSeuratQC_r(input, output, unique_id, adata_path=adata_path, assay=assay, min_genes=min_genes, max_genes=max_genes, min_UMI_count=min_UMI_count, max_UMI_count=max_UMI_count, percent_mt_max=percent_mt_max, percent_rb_min=percent_rb_min, resolution=resolution, dims=ro.r.seq(1, dims), doublet_rate=doublet_rate, regress_cell_cycle=ro.vectors.BoolVector([regress_cell_cycle])))
        if results[0] != ro.rinterface.NULL:
            default_assay = list(results[0])[0]
            assay_names = list(results[1])
            output = list(results[2])[0]
            adata_path = list(results[3])[0]
            ddl_assay_names = convert_from_r(list(results[4])[0])
        else:
            raise RuntimeError(f"Seurat QC failed.")

        if os.path.exists(adata_path):
            redislogger.info(unique_id, "Adding 3D UMAP to AnnData object.")
            adata = load_anndata(adata_path)
            sc.pp.neighbors(adata, n_neighbors=dims, n_pcs=n_pcs, random_state=0)
            adata_3D = sc.tl.umap(adata, random_state=0, 
                            init_pos="spectral", n_components=3, 
                            copy=True, maxiter=None)
            adata.obsm["X_umap_3D"] = adata_3D.obsm["X_umap"]
            adata.write_h5ad(adata_path, compression='gzip')
            adata_3D = None
        else:
            raise ValueError("AnnData file does not exist.")

    except Exception as e:
        redislogger.error(unique_id, f"An error happened while running Seurat QC: {e}.")

    return default_assay, assay_names, output, adata_path, adata, ddl_assay_names