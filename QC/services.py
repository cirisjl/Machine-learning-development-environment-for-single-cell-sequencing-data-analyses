# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import scanpy as sc
import dropkick as dk
from celery import Celery
# from fastapi.responses import FileResponse
import shutil
import hdf5plugin
from pydantic import BaseModel
from pathlib import Path
import numpy as np
import pandas as pd
import os
import os, time, sys
import sklearn
import scanpy as sc
import zipfile
from datetime import datetime
import xlrd

from datetime import datetime

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

import zipfile
import os.path

celeryApp = Celery('jobQueue', broker = 'pyamqp://')

def clearTemp(method):
    path = r"/tmp/{}/".format(method)
    now = time.time()

    for f in os.listdir(path):
        f = os.path.join(path, f)
        if os.stat(f).st_mtime < now - 60 * 60:
            # os.remove(os.path.join(path, f))
            shutil.rmtree(f)
    print('Removed temp files older than 1 hour.')

def read_data(method, filePath):
    clearTemp(method)
    ts = str(datetime.now().isoformat())

    tempFolder = '/tmp/{}/{}_'.format(method, method) + ts
    extension = filePath.split('.')[-1]
    if extension == 'zip': 
        with zipfile.ZipFile(filePath, 'r') as zip_ref:
            zip_ref.extractall(tempFolder)
        adata = sc.read_10x_mtx(
        tempFolder,  # the directory with the `.mtx` file
        var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
        cache=False)                              # write a cache file for faster subsequent reading

    else:
        if extension == 'csv':
            try:
                adata = sc.read_csv(filePath)
            except:
                adata = sc.read_text(filePath, delimiter='\t')
        elif 'h5' in extension:
            adata = sc.read_h5ad(filePath)
        elif 'xls' in extension:
            sheet_name = xlrd.open_workbook(filePath, on_demand=True).sheet_names()[0]
            adata = sc.read_excel(filePath,sheet_name)
        else:
            adata = sc.read_text(filePath, delimiter='\t')
        os.mkdir(tempFolder)
    
    return adata, tempFolder, ts


@celeryApp.task(name = 'ScanpyQC')
def processScanpyQC(filePath: str, min_genes: int, min_cells: int, n_genes_by_counts:int, pct_counts_mt:int):
    # clearTemp()
    # ts = str(datetime.now().isoformat())

    # tempFolder = '/tmp/scanpy/scanpy_' + ts
    # extension = filePath.split('.')[-1]
    # if extension == 'zip': 
    #     with zipfile.ZipFile(filePath, 'r') as zip_ref:
    #         zip_ref.extractall(tempFolder)
    #     adata = sc.read_10x_mtx(
    #     tempFolder,  # the directory with the `.mtx` file
    #     var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    #     cache=False)                              # write a cache file for faster subsequent reading

    # else:
    #     if extension == 'csv':
    #         try:
    #             adata = sc.read_csv(filePath)
    #         except:
    #             adata = sc.read_text(filePath, delimiter='\t')
    #     elif 'h5' in extension:
    #         adata = sc.read_h5ad(filePath)
    #     elif 'xls' in extension:
    #         sheet_name = xlrd.open_workbook(filePath, on_demand=True).sheet_names()[0]
    #         adata = sc.read_excel(filePath,sheet_name)
    #     else:
    #         adata = sc.read_text(filePath, delimiter='\t')
    #     os.mkdir(tempFolder)

    adata, tempFolder, ts = read_data('scanpy', filePath)

    adata.var_names_make_unique() 
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    adata = adata[adata.obs.n_genes_by_counts < n_genes_by_counts, :]
    adata = adata[adata.obs.pct_counts_mt < pct_counts_mt, :]
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata = adata[:, adata.var.highly_variable]
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    sc.tl.leiden(adata)
    sc.tl.leiden(adata, resolution=1.5, key_added="cluster2")

    sc.pl.umap(adata, color=['leiden','cluster2'], save = ts + '.png', show = False)

    os.mkdir(tempFolder + '/results')
    Path("figures/umap" + ts + '.png').rename(tempFolder + '/results/umap.png')

    adata.write_h5ad(
    tempFolder + '/results/processedAnnData.h5ad',
    compression=hdf5plugin.FILTERS["zstd"]
    )

    # create a ZipFile object
    with zipfile.ZipFile(tempFolder + '/output_' + ts + '.zip', 'w') as zipObj:
    # Iterate over all the files in directory
        for folderName, subfolders, filenames in os.walk(tempFolder + '/results'):
            for filename in filenames:
                #create complete filepath of file in directory
                filePath = os.path.join(folderName, filename)
                # Add file to zip
                zipObj.write(filePath, os.path.basename(filePath))

    # shutil.rmtree('/mnt/c/Users/Reddy/PycharmProjects/scanpy/scanpy_' + ts)
    # , background=os.remove('output_' + ts + '.zip')
    # return FileResponse(tempFolder + '/output_' + ts + '.zip')

@celeryApp.task(name = 'processDropkick')
def processDropkickQC(filePath: str, min_genes: int, min_cells: int, n_genes_by_counts:int, pct_counts_mt:int):
        
    # clearTemp()
    # ts = str(datetime.now().isoformat())

    # tempFolder = '/tmp/dropkick/dropkick_' + ts
    # extension = filePath.split('.')[-1]
    # if extension == 'zip': 
    #     with zipfile.ZipFile(filePath, 'r') as zip_ref:
    #         zip_ref.extractall(tempFolder)
    #     adata = sc.read_10x_mtx(
    #     tempFolder,  # the directory with the `.mtx` file
    #     var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    #     cache=False)                              # write a cache file for faster subsequent reading

    # else:
    #     if extension == 'csv':
    #         try:
    #             adata = sc.read_csv(filePath)
    #         except:
    #             adata = sc.read_text(filePath, delimiter='\t')
    #     elif 'h5' in extension:
    #         adata = sc.read_h5ad(filePath)
    #     elif 'xls' in extension:
    #         sheet_name = xlrd.open_workbook(filePath, on_demand=True).sheet_names()[0]
    #         adata = sc.read_excel(filePath,sheet_name)
    #     else:
    #         adata = sc.read_text(filePath, delimiter='\t')
    #     os.mkdir(tempFolder)

    adata, tempFolder, ts = read_data('dropkick', filePath)

    # adata = sc.read_10x_mtx(
    #     'data/tested/filtered_gene_bc_matrices/hg19',  # the directory with the `.mtx` file
    #     var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    #     cache=False)                              # write a cache file for faster subsequent reading

    adata = dk.recipe_dropkick(adata, n_hvgs=None, X_final="raw_counts")

    qc_plt = dk.qc_summary(adata, save_to = 'qcs_' + ts + '.png')

    # Run dropkick pipeline function

    adata_model = dk.dropkick(adata, n_jobs=5)
    adata_filtered = adata[adata.obs.dropkick_label == 'False'].copy()

    # Compare Dropkick output with Scanpy workflow
    score_plt = dk.score_plot(adata, save_to = 'score_' + ts + '.png')

    # # scanpy downstream processing on raw object
    # sc.pp.normalize_total(adata, target_sum=1e4)
    # sc.pp.log1p(adata)
    # sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    # adata = adata[:, adata.var.highly_variable]
    # sc.pp.regress_out(adata, ['total_counts'])
    # sc.pp.scale(adata, max_value=10)
    # sc.tl.pca(adata, svd_solver='arpack')
    # sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    # sc.tl.umap(adata)
    # sc.tl.leiden(adata)
    # sc.tl.leiden(adata, resolution=1.5, key_added="cluster2")

    # # Dropkick downstream processing
    # adata_filtered = dk.recipe_dropkick(adata_filtered, X_final="arcsinh_norm", filter=True, n_hvgs=2000, verbose=True)
    # sc.pp.neighbors(adata_filtered, n_neighbors=30, random_state=1, n_pcs=10)
    # sc.tl.leiden(adata_filtered)
    # sc.tl.umap(adata_filtered, random_state=1)

    # sc.pl.umap(adata, color=['dropkick_label'])

    # adata_filtered = adata[adata.obs.dropkick_label == 'False'].copy()
    # adata_filtered

    os.mkdir(tempFolder + '/results')
    Path("qcs_" + ts + '.png').rename(tempFolder + '/results/QC_Summary.png')
    Path("score_" + ts + '.png').rename(tempFolder + '/results/Dropkick_Score.png')
    
    return 
