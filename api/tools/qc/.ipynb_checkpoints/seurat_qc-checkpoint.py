import os
import numpy as np
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

# Defining the R script and loading the instance in Python
ro.r['source'](os.path.abspath(os.path.join(os.path.dirname(__file__), 'seurat_qc.R')))

def run_seurat_qc(input, output, save_anndata=True, assay='RNA', nFeature_min=200, nFeature_max=0, percent_mt_max=5, percent_rb_min=0, path_of_scrublet_calls=os.path.abspath(os.path.join(os.path.dirname(__file__), 'scrublet_calls.tsv')), dims=10, regress_cell_cycle=False):
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

    try:
        results = list(RunSeuratQC_r(input, output, save_anndata=ro.vectors.BoolVector([save_anndata]), assay=assay, nFeature_min=nFeature_min, nFeature_max=nFeature_max, percent_mt_max=percent_mt_max, percent_rb_min=percent_rb_min, path_of_scrublet_calls=path_of_scrublet_calls, dims=ro.r.seq(1, dims), regress_cell_cycle=ro.vectors.BoolVector([regress_cell_cycle])))
        default_assay = list(results[0])[0]
        assay_names = list(results[1])
        nCells = list(results[3])[0]
        nGenes = list(results[4])[0]
        genes  = list(results[5])
        cells = list(results[6])
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

    return default_assay, assay_names, metadata, nCells, nGenes, genes, cells, HVGsID, pca, tsne, umap