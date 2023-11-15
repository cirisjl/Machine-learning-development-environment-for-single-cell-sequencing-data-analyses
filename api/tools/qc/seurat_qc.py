import os
import numpy as np
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from tools.formating.formating import change_file_extension, load_anndata

# Ensure that pandas2ri is activated for automatic conversion
pandas2ri.activate()

# Defining the R script and loading the instance in Python
ro.r['source'](os.path.abspath(os.path.join(os.path.dirname(__file__), 'seurat_qc.R')))

def run_seurat_qc(input, assay='RNA', min_genes=200, max_genes=0, min_UMI_count=0, max_UMI_count=0, percent_mt_max=5, percent_rb_min=0, dims=10, regress_cell_cycle=False):
    RunSeuratQC_r = ro.globalenv['RunSeuratQC']
    
    default_assay = None
    assay_names = None

    if assay is None:
        assay = 'RNA'

    output = add_qc_result_suffix(input, assay)
    adata_path = change_file_extension(output, 'h5ad')

    print("paths")
    print(output)
    print(adata_path)

    try:
        results = list(RunSeuratQC_r(input, output, adata_path=adata_path, assay=assay, min_genes=min_genes, max_genes=max_genes, min_UMI_count=min_UMI_count, max_UMI_count=max_UMI_count,  percent_mt_max=percent_mt_max, percent_rb_min=percent_rb_min, dims=ro.r.seq(1, dims), regress_cell_cycle=ro.vectors.BoolVector([regress_cell_cycle])))
        if results[0] != ro.rinterface.NULL:
            default_assay = list(results[0])[0]
        assay_names = list(results[1])
        adata_path = list(results[2])[0]
        print(adata_path)

        if adata_path is not None:
            adata = load_anndata(adata_path)
            print("AnnData loaded successfully.")
        else:
            print("No AnnData path provided.")

    except Exception as e:
        print(e)

    return default_assay, assay_names, adata_path, adata


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