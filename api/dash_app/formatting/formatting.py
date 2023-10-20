import pandas as pd
import subprocess
import csv
import os
import scanpy as sc
from detect_delimiter import detect


# Convert Seurat/Single-Cell Experiment object to Anndata object and return the path of Anndata object
def convert_seurat_sce_to_anndata(path, assay='RNA'):
    import rpy2.robjects as ro

    # Get the absolute path of the current file
    current_file = os.path.abspath(__file__)

    # Construct the relative path to the desired file
    relative_path = os.path.join(os.path.dirname(current_file), 'formatting.R')

    # Get the absolute path of the desired file
    r_path = os.path.abspath(relative_path)

    with open(r_path, 'r') as r_source_file:
        r_source = r_source_file.read()

    # Evaluate the R script in the R environment
    ro.r(r_source)

    # Access the loaded R functions
    convert_seurat_sce_to_anndata = ro.globalenv['convert_seurat_sce_to_anndata']

    assay_names = None
    adata_path = None

    if path.endswith(".h5Seurat") or path.endswith(".h5seurat") or path.endswith(".rds"):
        try:
            results = convert_seurat_sce_to_anndata(path, assay=assay)
            adata_path = str(results.rx2('anndata_path'))
            assay_names = list(results[1])
        except Exception as e:
            print("Object format conversion is failed")
            print(e)

    return adata_path, assay_names