import os
import numpy as np
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

# Defining the R script and loading the instance in Python
ro.r['source'](os.path.abspath(os.path.join(os.path.dirname(__file__), 'seurat_qc.R')))

def run_seurat_qc(input, output, save_anndata=True, assay='RNA', min_genes=200, max_genes=0, min_UMI_count=0, max_UMI_count=0, percent_mt_max=5, percent_rb_min=0, path_of_scrublet_calls=os.path.abspath(os.path.join(os.path.dirname(__file__), 'scrublet_calls.tsv')), dims=10, regress_cell_cycle=False):
    RunSeuratQC_r = ro.globalenv['RunSeuratQC']
    default_assay = None
    assay_names = None
    metadata = None
    nCells = 0
    nGenes = 0
    genes = None
    cells = None
    HVGsID = None
    pca = None
    tsne = None
    umap = None
    adata_path = None

    try:
        results = list(RunSeuratQC_r(input, output, save_anndata=ro.vectors.BoolVector([save_anndata]), assay=assay, min_genes=min_genes, max_genes=max_genes, min_UMI_count=min_UMI_count, max_UMI_count=max_UMI_count,  percent_mt_max=percent_mt_max, percent_rb_min=percent_rb_min, path_of_scrublet_calls=path_of_scrublet_calls, dims=ro.r.seq(1, dims), regress_cell_cycle=ro.vectors.BoolVector([regress_cell_cycle])))
        if results[0] != ro.rinterface.NULL:
            default_assay = list(results[0])[0]
        assay_names = list(results[1])
        nCells = list(results[3])[0]
        nGenes = list(results[4])[0]
        if results[5] != ro.rinterface.NULL:
            genes  = list(results[5])
        if results[6] != ro.rinterface.NULL:
            cells = list(results[6])
        if results[11] != ro.rinterface.NULL:
            adata_path = list(results[11])[0]
        if results[7] != ro.rinterface.NULL:
            HVGsID = list(results[7])

        with localconverter(ro.default_converter + pandas2ri.converter):
            if results[2] != ro.rinterface.NULL: 
                metadata = ro.conversion.rpy2py(results[2])
            if results[8] != ro.rinterface.NULL:
                pca = ro.conversion.rpy2py(results[8])
            if results[9] != ro.rinterface.NULL:
                tsne = ro.conversion.rpy2py(results[9])
            if results[10] != ro.rinterface.NULL:
                umap = ro.conversion.rpy2py(results[10])
    except Exception as e:
        print(e)

    return default_assay, assay_names, metadata, nCells, nGenes, genes, cells, HVGsID, pca, tsne, umap, adata_path