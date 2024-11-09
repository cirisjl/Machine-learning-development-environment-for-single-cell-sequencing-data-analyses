import scanpy as sc
import os
import math
import hashlib
import sys
import subprocess
import numpy as np
import pandas as pd
from detect_delimiter import detect
import scipy.sparse as sp_sparse
from scipy.sparse import csr_matrix
from typing import Optional, Union
from string import ascii_letters
import csv
import gzip
import logging
import h5py
import jax
import jax.numpy as jnp
from collections import OrderedDict
from anndata import AnnData
from tools.visualization.plot import highest_expr_genes
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from tools.evaluation.clustering import clustering_scores
from tools.utils.gzip_str import *

from typing import Any, List, Optional
from attrdict import AttrDict
import json_numpy


# Ensure that pandas2ri is activated for automatic conversion
pandas2ri.activate()

# Defining the R script and loading the instance in Python
ro.r['source'](os.path.abspath(os.path.join(os.path.dirname(__file__), 'formating.R')))

try:
    from anndata._core.sparse_dataset import SparseDataset
except ImportError:
    # anndata >= 0.10.0
    from anndata._core.sparse_dataset import (
        BaseCompressedSparseDataset as SparseDataset,
    )


def load_anndata(path, annotation_path=None, dataset=None, assay='RNA', show_error=True, replace_invalid=False, isDashboard = False): # assay is optional and only for Seurat object
    # path = os.path.abspath(path)
    adata = None
    print(path)

    # if (os.path.isdir(path) and os.path.exists(os.path.join(path, "matrix.mtx")) and os.path.exists(
    #         os.path.join(path, "genes.tsv")) and os.path.exists(os.path.join(path, "barcodes.tsv"))) 
    #         or (os.path.isdir(path) and os.path.exists(os.path.join(path, "matrix.mtx.gz")) and os.path.exists(
    #         os.path.join(path, "genes.tsv.gz")) and os.path.exists(os.path.join(path, "barcodes.tsv.gz"))):
    if (os.path.isdir(path)):
        adata = sc.read_10x_mtx(path,
                             var_names='gene_symbols',  # use gene symbols for the variable names (variables-axis index)
                             cache=True)  # write a cache file for faster subsequent reading
    elif(os.path.exists(path)):
        # suffix = os.path.splitext(path)[-1]
        if path.endswith(".h5ad"):
            adata = sc.read_h5ad(path)
        elif path.endswith(".csv") or path.endswith(".tsv"):
            # print("Inside the loadAnndata CSV")
            print(detect_delimiter(path))
            # print("Inside the loadAnndata CSV 2")
            adata = sc.read_csv(path, delimiter=detect_delimiter(path))
            if annotation_path is not None:
                df_ann = pd.read_csv(annotation_path)
                df_ann = df_ann.set_index('cell', drop=False)
                adata.var = adata.var.join(df_ann)
            # print("Inside the loadAnndata CSV 3")
        elif path.endswith(".csv.gz") or path.endswith(".tsv.gz"):
            adata = sc.read_csv(path)
        elif path.endswith(".xlsx") or path.endswith(".xls"):
            adata = sc.read_excel(path, 0)
        # elif suffix == ".h5" and "pbmc" in path:
            # adata = sc.read_10x_h5(path)
        elif path.endswith(".h5"):
            try:
                adata = sc.read_10x_h5(path)
            except Exception as e:
                print(e)
                adata = sc.read_hdf(path, key=dataset)
        elif path.endswith(".loom"):
            adata = sc.read_loom(path)
        elif path.endswith(".mtx"):
            adata = sc.read_mtx(path)
        elif path.endswith(".txt") or path.endswith(".tab") or path.endswith(".data"):
            if replace_invalid:
                delimiter = detect_delimiter(path)
                adata = read_text_replace_invalid(path, delimiter)
            else:
                adata = sc.read_text(path, delimiter=detect_delim(path))
        elif path.endswith(".txt.gz"):
            if replace_invalid:
                adata = read_text_replace_invalid(path, "/t")
            else:
                adata = sc.read_text(path)
        elif path.endswith(".gz"):
            adata = sc.read_umi_tools(path)
        elif path.endswith(".h5Seurat") or path.endswith(".h5seurat") or path.endswith(".rds") or path.endswith(".Robj"):
            adata_path, assay_names, default_assay = convert_seurat_sce_to_anndata(path, assay=assay)
            if os.path.exists(adata_path):
                adata = sc.read_h5ad(adata_path)
        
    if adata is not None: 
        adata.obs = rename_col(adata.obs, 'n_counts')

    return adata


def change_file_extension(file_path, new_extension):
    # Check if the path is a directory
    if os.path.isdir(file_path):
        # If it's a directory, specify a default filename 'Anndata' with the new extension
        new_file_path = os.path.join(file_path, f"Anndata.{new_extension}")
    else:
        # If it's a file, proceed with changing the file's extension
        directory, base_filename = os.path.split(file_path)
        name, _ = os.path.splitext(base_filename)
        new_file_path = os.path.join(directory, f"{name}.{new_extension}").replace(" ", "_")

    return new_file_path


def convert_to_seurat_sce(input, output, format):
    import rpy2.rinterface_lib.callbacks as rcb
    import rpy2.robjects as ro
    ro.r['source'](os.path.abspath(os.path.join(os.path.dirname(__file__), 'formating.R')))

    if format == "Seurat":
        ConvertToSeurat_r = ro.globalenv['ConvertToSeurat']
        ConvertToSeurat_r(input, output)

    elif format == "SingleCellExperiment":
        ConvertToSCE_r = ro.globalenv['ConvertToSCE']
        ConvertToSCE_r(input, output)
    
    return output


def get_metadata_from_seurat(path):
    import rpy2.rinterface_lib.callbacks as rcb
    import rpy2.robjects as ro
    import anndata2ri
    from rpy2.robjects.packages import importr
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter

    rcb.logger.setLevel(logging.ERROR)
    ro.pandas2ri.activate()
    anndata2ri.activate()

    # Defining the R script and loading the instance in Python
    ro.r['source'](os.path.abspath(os.path.join(os.path.dirname(__file__), 'formating.R')))
    GetMetadataFromSeurat_r = ro.globalenv['GetMetadataFromSeurat']
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
    info = None

    try:
        results = list(GetMetadataFromSeurat_r(path))
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

        info = convert_from_r(results[11])

        # print(convert_from_r(results))
        # print(default_assay)
        # print(assay_names)

    except Exception as e:
        print("Error in get_metadata_from_seurat: ", e)

    return info, default_assay, assay_names, metadata, nCells, nGenes, genes, cells, HVGsID, pca, tsne, umap


def arr_to_list(arr:np.ndarray):
    arr_list = []
    for i in range(len(arr[0])):
        arr_list.append(arr[:, i].tolist())
    return arr_list


# def df_to_dict(df:pd.DataFrame):
#     df_dict = {}
#     df_columns = df.columns.values.tolist()
#     df_index = df.index.values.tolist()
#     df_dict['columns'] = df_columns
#     df_dict['index'] = df_index
    
#     for name in df_columns:
#         df_dict[name] = df[name].tolist()
        
#     return df_dict


def get_cell_metadata(adata, adata_path=None):
    cell_metadata = None
    obs_names = None
    nCells = 0
    nGenes = 0
    layers = None
    info = None
    adata_size = None
    cell_metadata_head = None
    embeddings = []

    if adata_path is not None and os.path.exists(adata_path):
        adata_size = file_size(adata_path)
    
    if adata is not None and isinstance(adata, AnnData):
        info = adata.__str__()
        layers = list(adata.layers.keys())
        obs = regularise_df(adata.obs)
        obs_names = obs.columns.values.tolist()
        nCells = adata.n_obs # Number of cells
        nGenes = adata.n_vars # Number of genes
        cell_metadata = df_to_dict(obs)
        cell_metadata_head = obs.dropna().head().to_dict()
        embedding_names = list(adata.obsm.keys()) # PCA, tSNE, UMAP
        for name in embedding_names:
            # embeddings.append({name: json_numpy.dumps(adata.obsm[name])})
            embeddings.append(name)
    
    return cell_metadata, cell_metadata_head, obs_names, nCells, nGenes, layers, info, adata_size, embeddings


def get_metadata_from_anndata(adata, pp_stage, process_id, process, method, parameters, md5, layer=None, adata_path=None, seurat_path=None, sce_path=None, cluster_label=None, scanpy_cluster='leiden'): 
    layers = None
    cell_metadata = None
    obs_names = None
    nCells = 0
    nGenes = 0
    genes = None
    gene_metadata = None
    embeddings = []
    umap = None
    tsne = None
    umap_3d = None
    # violin_plot = None
    # scatter_plot = None
    # highest_expr_genes_plot = None
    top_genes = None
    info = None
    pp_results = None
    adata_size = None
    seurat_size = None
    sce_size = None
    evaluation_results = []
    asw_score_leiden = None
    nmi_score_leiden = None
    ari_score_leiden = None
    asw_score_louvain = None
    nmi_score_louvain = None
    ari_score_louvain = None
    labels_pred_leiden = None
    labels_pred_louvain = None
    cluster_embedding = None
    description = f'{method} {process}' 

    if adata_path is not None and os.path.exists(adata_path):
        adata_size = file_size(adata_path)
    if seurat_path is not None and os.path.exists(seurat_path):
        seurat_size = file_size(seurat_path)
    if sce_path is not None and os.path.exists(sce_path):
        sce_size = file_size(sce_path)

    if adata is not None and isinstance(adata, AnnData):
        if layer is None:
            layer = "X"
            if cluster_label is not None and cluster_label in adata.obs.keys():
                cluster_label = adata.obs[cluster_label]
                if 'leiden' in adata.obs.keys() and layer+'_umap' in adata.obsm.keys():
                    labels_pred_leiden = adata.obs['leiden']
                if 'louvain' in adata.obs.keys() and layer+'_umap' in adata.obsm.keys():
                    labels_pred_louvain = adata.obs['louvain']
                cluster_embedding = adata.obsm[layer+'_umap']
        else:
            scanpy_cluster = layer + '_leiden'
            if cluster_label is not None and cluster_label in adata.obs.keys():
                cluster_label = adata.obs[cluster_label]
                if layer+'_leiden' in adata.obs.keys() and layer+'_umap' in adata.obsm.keys():
                    labels_pred_leiden = adata.obs[layer+'_leiden']
                if layer+'_louvain' in adata.obs.keys() and layer+'_umap' in adata.obsm.keys():
                    labels_pred_louvain = adata.obs[layer+'_louvain']
                cluster_embedding = adata.obsm[layer+'_umap']
        
        # if('cluster.ids' in adata.obs.keys()):
        #     scanpy_cluster = 'cluster.ids'
                
        info = adata.__str__()
        layers = list(adata.layers.keys())
        obs = regularise_df(adata.obs)
        obs_names = obs.columns.values.tolist()
        obs_dict = obs.to_dict('list') # Pandas dataframe
        obs_dict['index'] = obs.index.tolist()
        cell_metadata = gzip_dict(obs_dict)
        # cell_metadata = gzip_df(adata.obs) # pandas dataframe
        nCells = adata.n_obs # Number of cells
        nGenes = adata.n_vars # Number of genes
        if "highly_variable" in adata.var.keys():
            genes = gzip_list(adata.var[adata.var['highly_variable']==True].index.tolist())
            # gene_metadata = adata.var[adata.var['highly_variable']==True] # pandas dataframe
            var_dict = adata.var[adata.var['highly_variable']==True].to_dict('list') # Pandas dataframe
            var_dict['index'] = adata.var[adata.var['highly_variable']==True].index.tolist()
            gene_metadata = gzip_dict(var_dict)
        elif "vst.variable" in adata.var.keys():
            genes = gzip_list(adata.var[adata.var['vst.variable']==True].index.tolist())
            # gene_metadata = adata.var[adata.var['vst.variable']==True] # pandas dataframe
            var_dict = adata.var[adata.var['vst.variable']==True].to_dict('list') # Pandas dataframe
            var_dict['index'] = adata.var[adata.var['vst.variable']==True].index.tolist()
            gene_metadata = gzip_dict(var_dict)
        else:
            genes = gzip_list(adata.var_names.to_list()) # Gene IDs
            # gene_metadata = adata.var # pandas dataframe
            var_dict = adata.var.to_dict('list') # Pandas dataframe
            var_dict['index'] = adata.var.index.tolist()
            gene_metadata = gzip_dict(var_dict)

        embedding_names = list(adata.obsm.keys()) # PCA, tSNE, UMAP
        for name in embedding_names:
            # embeddings.append({name: json_numpy.dumps(adata.obsm[name])})
            embeddings.append(name)
        
        if layer != 'Pearson_residuals': # Normalize Pearson_residuals may create NaN values, which could not work with PCA
            if layer+'_umap' in adata.obsm.keys() and scanpy_cluster in adata.obs.keys():
                umap = json_numpy.dumps(adata.obsm[layer+'_umap'])
                # umap_plot = plot_UMAP(adata, layer=layer, clustering_plot_type=scanpy_cluster)
            elif layer+'_umap' in adata.obsm.keys():
                umap = json_numpy.dumps(adata.obsm[layer+'_umap'])
                # umap_plot = plot_UMAP(adata, layer=layer)
            
            if layer+'_umap_3D' in adata.obsm.keys() and scanpy_cluster in adata.obs.keys():
                umap_3d = json_numpy.dumps(adata.obsm[layer+'_umap_3D'])
                # umap_plot_3d = plot_UMAP(adata, layer=layer, clustering_plot_type=scanpy_cluster, n_dim=3)
            elif layer+'_umap_3D' in adata.obsm.keys():
                umap_3d = json_numpy.dumps(adata.obsm[layer+'_umap_3D'])
                # umap_plot_3d = plot_UMAP(adata, layer=layer, n_dim=3)

            if layer+'_tsne' in adata.obsm.keys():
                tsne = json_numpy.dumps(adata.obsm[layer+'_tsne'])

        if process == 'QC':
            # violin_plot = gzip_str(plot_violin(adata))
            # scatter_plot = gzip_str(plot_scatter(adata))
            # violin_plot = plot_violin(adata)
            # scatter_plot = plot_scatter(adata)
            if nCells < 12000: # If the dataset is too large, then skip the highest expressed genes plot
                counts_top_genes, columns = highest_expr_genes(adata)
                top_genes = {"counts_top_genes": json_numpy.dumps(counts_top_genes), "columns": columns}

        if cluster_label is not None:
            if labels_pred_leiden is not None:
                asw_score_leiden, nmi_score_leiden, ari_score_leiden = clustering_scores(cluster_label, labels_pred_leiden, cluster_embedding)
                evaluation_results.append(
                    {
                        "leiden": {
                            "asw_score": asw_score_leiden,
                            "nmi_score": nmi_score_leiden,
                            "ari_score": ari_score_leiden
                        }
                    }
                )
            if labels_pred_louvain is not None:
                asw_score_louvain, nmi_score_louvain, ari_score_louvain = clustering_scores(cluster_label, labels_pred_louvain, cluster_embedding)
                evaluation_results.append(
                    {
                        "louvain": {
                            "asw_score": asw_score_louvain,
                            "nmi_score": nmi_score_louvain,
                            "ari_score": ari_score_louvain
                        }
                    }
                )
            evaluation_results = {
                "leiden": {
                    "asw_score": asw_score_leiden,
                    "nmi_score": nmi_score_leiden,
                    "ari_score": ari_score_leiden
                },
                "louvain": {
                    "asw_score": asw_score_louvain,
                    "nmi_score": nmi_score_louvain,
                    "ari_score": ari_score_louvain
                }
            }

        pp_results = {
            "process_id": process_id,
            "description": description,
            "md5": md5,
            "stage": pp_stage,
            "process": process,
            "method": method,
            "parameters": parameters,
            "info": info,
            "adata_path": adata_path,
            "seurat_path": seurat_path,
            "sce_path": sce_path,
            "adata_size": adata_size,
            "seurat_size": seurat_size,
            "sce_size": sce_size,
            "layer": layer,
            "layers": layers,
            "obs_names": obs_names,
            "cell_metadata": cell_metadata,
            "gene_metadata": gene_metadata,
            "nCells": nCells,
            "nGenes": nGenes,
            "genes": genes,
            "embeddings": embeddings,
            "umap": umap,
            "umap_3d": umap_3d,
            "tsne": tsne,
            "highest_expr_genes": top_genes,
            # "violin_plot": violin_plot,
            # "scatter_plot": scatter_plot,
            # "highest_expr_genes_plot": highest_expr_genes_plot,
            "evaluation_results": evaluation_results
            }
        
    return pp_results


def file_size(path): # MB
    return round(os.path.getsize(path)/(1024*1024), 2)


# Convert Seurat/Single-Cell Experiment object to Anndata object and return the path of Anndata object
def convert_seurat_sce_to_anndata(path, assay='RNA'):

    if assay is None:
        assay = 'RNA'

    # Access the loaded R functions
    ConvertSeuratSCEtoAnndata_r = ro.globalenv['ConvertSeuratSCEtoAnndata']

    assay_names = None
    adata_path = None
    default_assay = None

    if path.endswith(".h5Seurat") or path.endswith(".h5seurat") or path.endswith(".rds") or path.endswith(".Robj"):
        try:
            print("convert_seurat_sce_to_anndata")
            print(assay)
            results = list(ConvertSeuratSCEtoAnndata_r(path, assay=assay))
            if results[0] is not None and results[0] != ro.rinterface.NULL:
                default_assay = list(results[0])[0]
            else:
                default_assay = None  # or a sensible default like 'RNA'

            if results[1] is not None and results[1] != ro.rinterface.NULL:
                assay_names = list(results[1])
            else:
                assay_names = []

            if results[2] is not None and results[2] != ro.rinterface.NULL:
                adata_path = list(results[2])[0]
            else:
                adata_path = None     

        except Exception as e:
            print("Seurat/SCE to Anndata is failed")
            print(e)

    return adata_path, assay_names, default_assay


def anndata_to_csv(adata, output_path, layer=None, compress=False):
    counts = None

    if layer is None:
        if type(adata.X) != "numpy.ndarray":
            counts = adata.X.toarray()
        else:
            counts = adata.X
    else:
        if type(adata.layers[layer]) != "numpy.ndarray":
            counts = adata.layers[layer].toarray()
        else:
            counts = adata.layers[layer]

    pd.DataFrame(data=counts, index=adata.obs_names, columns=adata.var_names).to_csv(output_path)
    return output_path


def load_anndata_to_csv(input, csv_path, layer=None, show_error=True, dataset=None, compress=False):
    adata = None
    counts = None

    try:
        adata = load_anndata(input, dataset)
        print(adata)
        adata_path = input
    except Exception as e:
        print("File format is not supported.")
        if show_error: print(e)
        return None, None, None

    if layer is None:
        counts = adata.X
    elif layer in adata.layers.keys():
        counts = adata.layers[layer]       
    else:
        print("Layer is not found in AnnData object.")
        return None, None, None

    csv_path = anndata_to_csv(adata, csv_path, layer=layer, compress=compress)

    return adata, counts, csv_path


def detect_delim(path):
    # look at the first ten thousand bytes to guess the character encoding
    with open(path, 'rb') as file:
        rawdata = file.read(10000)
        rawdata = rawdata.decode('utf-8')
        delimiter = detect(rawdata, whitelist=[' ', ',', ';', ':', '|', '\t'])
        return delimiter


def detect_delimiter(file_path):
    with open(file_path, 'r') as file:
        # Read the first line of the file to detect the delimiter
        first_line = file.readline()
        dialect = csv.Sniffer().sniff(first_line)
        return dialect.delimiter


# def output_path_check(dataset, output, method = '', format = "AnnData"):
#     output = os.path.abspath(output)
#     if method != '': method = '_' + method
    
#     if not os.path.exists(os.path.dirname(output)):
#         os.makedirs(os.path.dirname(output))

#     if os.path.isdir(output) and format == "AnnData":
#         output = os.path.join(output, dataset + method + ".h5ad")
#         print("The output path is a directory, adding output file " + dataset + method + ".h5ad to the path.")
#     elif os.path.isdir(output) and format == "SingleCellExperiment":
#         output = os.path.join(output, dataset + method + ".rds")
#         print("The output path is a directory, adding output file " + dataset + method + ".rds to the path.")
#     elif os.path.isdir(output) and format == "Seurat":
#         output = os.path.join(output, dataset + method + ".h5Seurat")
#         print("The output path is a directory, adding output file " + dataset + method + ".h5Seurat to the path.")
#     elif os.path.isfile(output) and format == "AnnData" and os.path.splitext(output)[-1] != ".h5ad":
#         output.replace(os.path.splitext(output)[-1], method + ".h5ad")
#         print("The suffix is incorrect, changing it to '.h5ad'.")
    
#     return output


def get_output_path(path, process_id='', dataset=None, method='', format="AnnData", compress=False):
    output = os.path.abspath(path)
    method = '_' + method if method != '' else ''
    output_path = None
    directory = output
    base_name = os.path.basename(path)

    if not os.path.exists(output):
        os.makedirs(output)

    if os.path.isdir(output):
        if dataset is None:
            dataset = base_name
        if format == "AnnData":
            output_path = os.path.join(output, process_id, dataset + method + ".h5ad")
            print("The output path is a directory, adding output file " + dataset + method + ".h5ad to the path.")
        elif format == "SingleCellExperiment":
            output_path = os.path.join(output, process_id, dataset + method + ".rds")
            print("The output path is a directory, adding output file " + dataset + method + ".rds to the path.")
        elif format == "Seurat":
            output_path = os.path.join(output, process_id, dataset + method + ".h5seurat")
            print("The output path is a directory, adding output file " + dataset + method + ".h5seurat to the path.")
        elif format == "CSV":
            if compress == True:
                output_path = os.path.join(output, process_id, dataset + method + ".csv.gz")
                print("The output path is a directory, adding output file " + dataset + method + ".csv.gz to the path.")
            else:
                output_path = os.path.join(output, process_id, dataset + method + ".csv")
                print("The output path is a directory, adding output file " + dataset + method + ".csv to the path.")
    else:
        if format == "AnnData":
            output_path = os.path.join(directory, process_id, base_name.replace(os.path.splitext(output)[-1], method + ".h5ad"))
        elif format == "SingleCellExperiment":
            output_path = os.path.join(directory, process_id, base_name.replace(os.path.splitext(output)[-1], method + ".rds"))
        elif format == "Seurat":
            output_path = os.path.join(directory, process_id, base_name.replace(os.path.splitext(output)[-1], method + ".h5seurat"))
        elif format == "CSV":
            if compress == True:
                output_path = os.path.join(directory, process_id, base_name.replace(os.path.splitext(output)[-1], method + ".csv.gz"))
            else:
                output_path = os.path.join(directory, process_id, base_name.replace(os.path.splitext(output)[-1], method + ".csv"))

    if not os.path.exists(os.path.dirname(output_path)):
        os.makedirs(os.path.dirname(output_path))
    
    output_path = output_path.replace(" ", "_")
    print(output_path)

    return output_path


def get_report_path(dataset, output, method):
    output = os.path.abspath(output)
    method = '_' + method if method else ''
    report_path = None

    if not os.path.exists(os.path.dirname(output)):
        os.makedirs(os.path.dirname(output))

    report_path = output.replace(os.path.splitext(output)[-1], method + "_report.html")
    
    return report_path


def list_py_to_r(list):
    list = [x.upper() for x in list if isinstance(x,str)]
    return 'c(' + ','.join(list) + ')'


def methods_list(list):
    list = [x.upper() for x in list if isinstance(x,str)]
    return ','.join(list)


def list_to_string(list):
    list = [x.upper() for x in list if isinstance(x, str)]
    return ','.join(list)


def list_to_string_default(list):
    list = [x for x in list if isinstance(x, str)]
    return ','.join(list)


def convert_gz_to_txt(gz_file_path, txt_file_path):
  with gzip.open(gz_file_path, 'rb') as f_in:
    with open(txt_file_path, 'w') as f_out:
      f_out.write(f_in.read().decode())


def read_text_replace_invalid(file_path, delimiter):
    if file_path.endswith(".gz"):
        file_name_without_extension = os.path.splitext(os.path.basename(file_path))[0]
        # Create the new file path with a .txt extension
        new_file_path = os.path.join(os.path.dirname(file_path), file_name_without_extension + '.txt')
        convert_gz_to_txt(file_path, new_file_path)
        df = pd.read_csv(new_file_path, sep="\t", on_bad_lines='skip', index_col=0)
    else:
        df = pd.read_csv(file_path, delimiter=delimiter, on_bad_lines='skip', index_col=0)
    
    df = df.apply(pd.to_numeric, errors='coerce')
    return sc.AnnData(df)


def read_text(file_path):
    if file_path.endswith(".gz"):
        
        file_name_without_extension = os.path.splitext(os.path.basename(file_path))[0]

        # Create the new file path with a .txt extension
        new_file_path = os.path.join(os.path.dirname(file_path), file_name_without_extension + '.txt')
        convert_gz_to_txt(file_path, new_file_path)

        df = pd.read_csv(new_file_path, sep="\t", on_bad_lines='skip', index_col=0)
        return sc.AnnData(df)
    else:
        delimiter = detect_delimiter(file_path)
        df = pd.read_csv(file_path, delimiter=delimiter, on_bad_lines='skip', index_col=0)
        return sc.AnnData(df)
    

def load_invalid_adata(file_path, replace_nan):
    if file_path.endswith(".gz"):
        file_name_without_extension = os.path.splitext(os.path.basename(file_path))[0]

        # Create the new file path with a .txt extension
        new_file_path = os.path.join(os.path.dirname(file_path), file_name_without_extension + '.txt')
        convert_gz_to_txt(file_path, new_file_path)
        df = pd.read_csv(new_file_path, sep="\t", on_bad_lines='skip', index_col=0)
    else:
        delimiter = detect_delimiter(file_path)
        df = pd.read_csv(file_path, delimiter=delimiter, on_bad_lines='skip', index_col=0)

    invalid_rows = df.apply(pd.to_numeric, errors='coerce').isnull().any(axis=1)
    invalid_columns = df.columns[df.apply(pd.to_numeric, errors='coerce').isnull().any()]
    invalid_df = df.loc[invalid_rows, invalid_columns]

    return sc.AnnData(invalid_df)


def is_normalized(expression_matrix, min_genes):
    if (isinstance(expression_matrix, csr_matrix)):
        expression_matrix = expression_matrix.toarray()

    if (isinstance(expression_matrix, np.ndarray)):
        if np.min(expression_matrix) < 0 or np.max(expression_matrix) < min_genes:
            return True
        else:
            return False
        

def check_nonnegative_integers(
    data: Union[pd.DataFrame, np.ndarray, sp_sparse.spmatrix, h5py.Dataset],
    n_to_check: int = 20,
):
    """Approximately checks values of data to ensure it is count data."""
    # for backed anndata
    if isinstance(data, h5py.Dataset) or isinstance(data, SparseDataset):
        data = data[:100]

    if isinstance(data, np.ndarray):
        data = data
    elif issubclass(type(data), sp_sparse.spmatrix):
        data = data.data
    elif isinstance(data, pd.DataFrame):
        data = data.to_numpy()
    else:
        raise TypeError("data type not understood")

    ret = True
    if len(data) != 0:
        inds = np.random.choice(len(data), size=(n_to_check,))
        check = jax.device_put(data.flat[inds], device=jax.devices("cpu")[0])
        negative, non_integer = _is_not_count_val(check)
        ret = not (negative or non_integer)
    return ret


@jax.jit
def _is_not_count_val(data: jnp.ndarray):
    negative = jnp.any(data < 0)
    non_integer = jnp.any(data % 1 != 0)

    return negative, non_integer


def get_file_md5(path: str, split_num=256, get_byte=8):

    if not isinstance(split_num, int) or split_num <= 0:
        raise TypeError("split_num must be a positive none-zero integer!")
    if not isinstance(get_byte, int) or get_byte <= 0:
        raise TypeError("get_byte must be a positive none-zero integer!")
    if not os.path.exists(path):
        raise TypeError("%s does not exist!" % path)
    if os.path.isdir(path):
        raise TypeError("%s is a folder, while path should be a file!" % path)
    
    size = round(os.path.getsize(path), 2)
    # For a small file (equal to or less than 2M), caculate the MD5 values directly.
    if size < split_num * get_byte:
        # Read the file
        with open(path, 'rb') as f1:
            f1 = f1.read()
        cipher = hashlib.md5()
        cipher.update(str(split_num).encode('utf-8'))
        cipher.update(f1)
        cipher.update(str(get_byte).encode('utf-8'))
        return cipher.hexdigest()
    # For a large file, split the file in to several segment, and then sum the MD5 value. 
    mean_size = size // split_num
    cipher = hashlib.md5()
    # Position
    place = 0
    with open(path, 'rb') as f1:
        for i in range(split_num):
            f1.seek(place)
            res = f1.read(get_byte)
            cipher.update(res)
            place = place + mean_size

    return cipher.hexdigest()


def get_md5(path:str):
    md5 = []
    if not os.path.exists(path):
        raise TypeError("%s does not exist!" % path)
    if os.path.isdir(path):
        for file in os.listdir(path):
            md5.append(get_file_md5(path+file))
    else:
        md5.append(get_file_md5(path))
    return md5


def convert_df_dates_from_r(df: pd.DataFrame, date_cols: Optional[List[str]] = None) -> pd.DataFrame:
    """ convert given date columns into pandas datetime with UTC timezone
    Args:
        df (pd.DataFrame): The pandas datframe
        date_cols (list[str], optional): _description_. Defaults to None.
    Returns:
        pd.DataFrame: The dataframe with the converted
    """
    result = df.copy()
    if date_cols is not None:
        for col in (set(date_cols) & set(result.columns)):
            result[col] = pd.to_datetime(
                result[col], unit='D', origin='1970-1-1').dt.tz_localize('UTC')
    return result


def convert_to_r(item: Any) -> Any:
    """ cpnverts python object into rpy2 format
    Args:
        item (Any): native python object
    Returns:
        Any: rpy2 object
    """
    if item is None:
        return ro.r("NULL")
    elif isinstance(item, pd.DataFrame):
        with localconverter(ro.default_converter + pandas2ri.converter):
            result = ro.conversion.py2rpy(item)
        return result
    elif isinstance(item, np.ndarray):
        return ro.FloatVector(item)
    elif isinstance(item, (AttrDict, pd.Series)):
        return convert_to_r(dict(item))
    elif isinstance(item, dict):
        temp = {k: convert_to_r(v) for k, v in item.items()
                if v is not None}
        temp = {k: v for k, v in temp.items() if v is not None}
        return ro.ListVector(temp)
    elif isinstance(item, set):
        return convert_to_r(list(item))
    elif isinstance(item, (list, tuple, pd.Index)):
        if len(item) == 0:
            return None
        if isinstance(item[0], float):
            return ro.FloatVector(item)
        if isinstance(item[0], (int, np.int0)):
            return ro.IntVector(item)
        return ro.StrVector([str(i) for i in item])
    else:
        return item


def convert_from_r(item: Any, date_cols: Optional[List[str]] = None, name: str = '', reserve_plots: bool = True) -> Any:
    """convert rpy object into python native object
    Args:
        item (Any): rpy2 object to convert
        date_cols (list[str], optional): define the date colums in R dataframe , in order to conevrt them into pandas datetime. Defaults to None.
        name (str, optional): name of the object to convert, (not required for external use). Defaults to ''.
        reserve_plots (bool, optional): if True prserve rpy2 ListVector as rpy2 ListVector if name conatains plot,
                                        in order to be able to ouput ggplot plots. Defaults to True.
    Returns:
        Any: the converted item
    """
    result = item
    remove_list: bool = True
    if item == ro.vectors.NULL:
        return None
    elif 'plot' in name and isinstance(item, ro.vectors.ListVector) and reserve_plots:
        return item
    elif isinstance(item, (ro.environments.Environment,
                           ro.Formula)):
        return None
    elif isinstance(item, ro.vectors.DataFrame):
        with localconverter(ro.default_converter + pandas2ri.converter):
            result = ro.conversion.rpy2py(item)
        result = convert_df_dates_from_r(result, date_cols)
        remove_list = False
    elif isinstance(item, (ro.vectors.StrVector,
                           ro.vectors.FloatVector,
                           ro.vectors.BoolVector,
                           ro.vectors.IntVector)):
        result = tuple(item)
    elif isinstance(item, ro.vectors.ListVector):
        if item.names == ro.vectors.NULL:
            return None
        result = {}
        remove_list = False
        if len(item) > 0:
            result = dict(zip(item.names, list(item)))
            for k, v in result.items():
                result[k] = convert_from_r(v, date_cols, name=k)
    if '__len__' in result.__dir__() and len(result) == 1 and remove_list:
        result = result[0]
    return result


# def load_annData_dash(path, replace_invalid=False):
#     show_error=True
#     dataset = None
#     # path = os.path.abspath(path)
#     adata = None
    
#     if (os.path.isdir(path)):
#         adata = sc.read_10x_mtx(path,
#                              var_names='gene_symbols',  # use gene symbols for the variable names (variables-axis index)
#                              cache=True)  # write a cache file for faster subsequent reading
#     elif(os.path.exists(path)):
#         # suffix = os.path.splitext(path)[-1]
#         if path.endswith(".h5ad"):
#             adata = sc.read_h5ad(path)
#         elif path.endswith(".csv") or path.endswith(".tsv"):
#             # print("Inside the loadAnndata CSV")
#             print(detect_delimiter(path))
#             # print("Inside the loadAnndata CSV 2")
#             adata = sc.read_csv(path, delimiter=detect_delimiter(path))
#             # print("Inside the loadAnndata CSV 3")
#         elif path.endswith(".csv.gz") or path.endswith(".tsv.gz"):
#             data = sc.read_csv(path)
#         elif path.endswith(".xlsx") or path.endswith(".xls"):
#             adata = sc.read_excel(path, 0)
#         # elif suffix == ".h5" and "pbmc" in path:
#             # adata = sc.read_10x_h5(path)
#         elif path.endswith(".h5"):
#             try:
#                 adata = sc.read_10x_h5(path)
#             except Exception as e:
#                 print(e)
#                 adata = sc.read_hdf(path, key=dataset)
#         elif path.endswith(".loom"):
#             adata = sc.read_loom(path)
#         elif path.endswith(".mtx"):
#             adata = sc.read_mtx(path)
#         elif path.endswith(".txt") or path.endswith(".tab") or path.endswith(".data"):
#             delimiter = detect_delimiter(path)
#             if replace_invalid:
#                 adata = read_text_replace_invalid(path, delimiter)
#                 print(adata)
#                 print(adata.var_names[:10])
#                 print(adata.obs_names[:10])
#             else:
#                 adata = sc.read_text(path, delimiter=detect_delimiter(path))      
#         elif path.endswith(".txt.gz"):
#             if replace_invalid:
#                 adata = read_text_replace_invalid(path, "/t")
#             else:
#                 adata = sc.read_text(path)      
#         elif path.endswith(".gz"):
#             adata = sc.read_umi_tools(path)
#         elif path.endswith(".h5Seurat") or path.endswith(".h5seurat") or path.endswith(".rds"):
#             try:
#                 current_file = os.path.abspath(__file__)
#                 # Construct the relative path to the desired file
#                 relative_path = os.path.join(os.path.dirname(current_file), 'convert_to_anndata.Rmd')

#                 # Get the absolute path of the desired file
#                 operation_path = os.path.abspath(relative_path)
#                 report_path = os.path.join(os.path.dirname(path), "file_conversion_report.html")
#                 adata_path = os.path.splitext(path)[0] + '.h5ad'
                
#                 if os.path.exists(adata_path):
#                     adata = sc.read_h5ad(adata_path)
#                 else:
#                     s = subprocess.call(["R -e \"rmarkdown::render('" + operation_path + "', params=list(path='" + str(path) + "'), output_file='" + report_path + "')\""], shell = True)
#                     print(s)
#                     adata = sc.read_h5ad(adata_path)

#             except Exception as e:
#                 print("Object format conversion is failed")
#                 if show_error: print(e)

#     return adata

# Remove NA and single value columns
def regularise_df(df):
    df = df.dropna(axis=1, how='all')
    res = df
    for col in df.columns:
        if len(df[col].unique()) == 1:
            res = res.drop(col,axis=1)
    return res


# Drop numerical columns for cell type annotation
def drop_num_col(df):
    col_to_drop = []

    for i, v in df.dtypes.items():
        if 'float' in str(v) or 'int' in str(v) or is_number(df[i][0]):
            col_to_drop.append(i)
    
    if len(col_to_drop) > 0:
        df = df.drop(columns=col_to_drop)
    
    return df


def rename_col(df, pattern):
    import re
    for name in df.columns.values:
        if re.match(pattern, name):
            df.rename(columns={name: pattern}, inplace=True)
            break
    return df


# Convert dataframe to dict for cell type selection
def df_to_dict(df):
    df_dict = {}
    # df = regularise_df(df)
    df = drop_num_col(df)
    for col in df.columns:
        col_dict = {}
        if len(df[col].unique()) < 500 and isinstance(df[col][0], str): # The known human cell types are less than 500
            col_dict[col] = [x for x in df[col].unique().tolist() if isinstance(x, str) or not math.isnan(x)]
            # col_dict[col] = df[col].unique().tolist()
            df_dict[col] = col_dict
    
    return df_dict


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        pass
 
    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass
 
    return False


def save_anndata(adata, output):
    if not isinstance(adata.X, csr_matrix):
        adata.X = csr_matrix(adata.X)
    if adata.raw is not None:
        if not isinstance(adata.raw.X, csr_matrix):
            adata.raw.X = csr_matrix(adata.raw.X)
    adata.write_h5ad(output, compression='gzip')
    
    return output