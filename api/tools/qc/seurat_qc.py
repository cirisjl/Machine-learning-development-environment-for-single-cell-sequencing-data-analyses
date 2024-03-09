import os
import scanpy as sc
import numpy as np
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from tools.formating.formating import change_file_extension, load_anndata, convert_from_r, convert_to_r

# Ensure that pandas2ri is activated for automatic conversion
pandas2ri.activate()

# Defining the R script and loading the instance in Python
ro.r['source'](os.path.abspath(os.path.join(os.path.dirname(__file__), 'seurat_qc.R')))

def run_seurat_qc(input, assay='RNA', min_genes=200, max_genes=0, min_UMI_count=0, max_UMI_count=0, percent_mt_max=5, percent_rb_min=0, resolution=0.5, dims=10, doublet_rate=0.075, regress_cell_cycle=False):
    RunSeuratQC_r = ro.globalenv['RunSeuratQC']
    adata = None
    default_assay = None
    assay_names = None
    ddl_assay_names = False

    if assay is None:
        assay = 'RNA'

    output = add_qc_result_suffix(input, assay)
    adata_path = change_file_extension(output, 'h5ad')

    try:
        results = list(RunSeuratQC_r(input, output, adata_path=adata_path, assay=assay, min_genes=min_genes, max_genes=max_genes, min_UMI_count=min_UMI_count, max_UMI_count=max_UMI_count, percent_mt_max=percent_mt_max, percent_rb_min=percent_rb_min, resolution=resolution, dims=ro.r.seq(1, dims), doublet_rate=doublet_rate, regress_cell_cycle=ro.vectors.BoolVector([regress_cell_cycle])))
        if results[0] != ro.rinterface.NULL:
            default_assay = list(results[0])[0]
            assay_names = list(results[1])
            adata_path = list(results[2])[0]
            ddl_assay_names = convert_from_r(list(results[3])[0])

        print(adata_path)

        if os.path.exists(adata_path):
            adata = load_anndata(adata_path)
            adata_3D = sc.tl.umap(adata, random_state=0, 
                            init_pos="spectral", n_components=3, 
                            copy=True, maxiter=None)
            adata.obsm["X_umap_3D"] = adata_3D.obsm["X_umap"]
            adata.write_h5ad(adata_path)
            adata_3D = None
        else:
            print("AnnData file does not exist.")

    except Exception as e:
        print(e)

    return default_assay, assay_names, adata_path, adata, output, ddl_assay_names


def add_qc_result_suffix(input_path, assay):
    # Split the input path into the directory, filename, and extension
    directory, filename = os.path.split(input_path)
    base_name, extension = os.path.splitext(filename)

    # If assay is 'RNA', add underscore before 'qc_result', otherwise, include assay directly
    if assay == 'RNA':
        new_base_name = base_name + "_qc_result"
    else:
        new_base_name = base_name + "_qc_result_" + assay

    # Reconstruct the output path
    output_path = os.path.join(directory, new_base_name + extension)
    return output_path