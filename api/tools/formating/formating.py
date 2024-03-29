import scanpy as sc
import os
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
from tools.visualization.plot import plot_UMAP, plot_scatter, plot_highest_expr_genes, plot_violin
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

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
        new_file_path = os.path.join(directory, f"{name}.{new_extension}")

    return new_file_path


def convert_to_seurat_sce(input, output, format):
    import rpy2.rinterface_lib.callbacks as rcb
    import rpy2.robjects as ro
    ro.r['source'](os.path.abspath(os.path.join(os.path.dirname(__file__), 'formating.R')))

    if format == "Seurat":
        ConvertToSeurat_r = ro.globalenv['ConvertToSeurat']
        output = ConvertToSeurat_r(input, output)
        output = convert_from_r(output)

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


def get_metadata_from_anndata(adata, pp_stage, process_id, process, method, parameters, layer=None, adata_path=None, seurat_path=None, sce_path=None):
    layers = None
    cell_metadata_obs = None
    nCells = 0
    nGenes = 0
    genes = None
    cells = None
    gene_metadata = None
    embeddings = []
    umap_plot = None
    umap_plot_3d = None
    violin_plot = None
    scatter_plot = None
    highest_expr_genes_plot = None
    info = None
    pp_results = None

    if adata is not None and isinstance(adata, AnnData):
        info = adata.__str__()
        layers = list(adata.layers.keys())
        cell_metadata_obs = adata.obs # pandas dataframe
        nCells = adata.n_obs # Number of cells
        nGenes = adata.n_vars # Number of genes
        genes = adata.var_names.to_list() # Cell IDs
        cells = adata.obs_names.to_list() # Gene IDs
        gene_metadata = adata.var # pandas dataframe
        embedding_names = list(adata.obsm.keys()) # PCA, tSNE, UMAP
        for name in embedding_names:
            embeddings.append({name: json_numpy.dumps(adata.obsm[name])})

        umap_plot = plot_UMAP(adata, layer=layer)
        umap_plot_3d = plot_UMAP(adata, layer=layer, n_dim=3)
        violin_plot = plot_violin(adata)
        scatter_plot = plot_scatter(adata)
        highest_expr_genes_plot = plot_highest_expr_genes(adata)

        pp_results = {
            "process_id": process_id,
            "stage": pp_stage,
            "process": process,
            "method": method,
            "parameters": parameters,
            "info": info,
            "adata_path": adata_path,
            "seurat_path": seurat_path,
            "sce_path": sce_path,
            "layers": layers,
            "cell_metadata_obs": cell_metadata_obs.to_dict(),
            "gene_metadata": gene_metadata.to_dict(),
            "nCells": nCells,
            "nGenes": nGenes,
            "genes": genes,
            "cells": cells,
            "embeddings": embeddings,
            "umap_plot": umap_plot,
            "umap_plot_3d": umap_plot_3d,
            "violin_plot": violin_plot,
            "scatter_plot": scatter_plot,
            "highest_expr_genes_plot": highest_expr_genes_plot
            }
        
    return pp_results


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


def anndata_to_csv(adata, output_path, layer = None):
    counts = None
    print("to CSV")
    print(adata)

    if layer is None:
        if type(adata.layers[layer]) != "numpy.ndarray":
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


def load_anndata_to_csv(input, output, layer=None, show_error=True, dataset=None):
    adata = None
    adata_path = None
    counts = None

    if os.path.exists(output):
        try:
            adata = load_anndata(output, dataset)
            adata_path = output
        except Exception as e:
            print("File format is not supported.")
            if show_error: print(e)
            return None, None, None
    else:
        try:
            print("Inside else , read from input path")
            print(input)
            adata = load_anndata(input, dataset)
            print(adata)
            adata_path = input
        except Exception as e:
            print("File format is not supported.")
            if show_error: print(e)
            return None, None, None

    if layer is None:
        counts = adata.X
        print("Layer is none")
        print(counts)
    elif layer in adata.layers.keys():
        counts = adata.layers[layer]       
    else:
        print("Layer is not found in AnnData object.")
        return None, None, None

    csv_path = adata_path.replace(".h5ad", ".csv")
    csv_path = anndata_to_csv(adata, csv_path, layer=layer)

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


def get_output_path(path, process_id='', dataset=None, method = '', format = "AnnData"):
    output = os.path.abspath(path)
    method = '_' + method if method else ''
    directory, base_name = os.path.split(output.rstrip('/'))
    output_path = None

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
            output_path = os.path.join(directory, process_id, base_name.replace(os.path.splitext(output)[-1], method + ".csv"))
    
    return output_path


def get_report_path(dataset, output, method):
    output = os.path.abspath(output)
    report_path = None
    if os.path.isdir(output):
        report_path = os.path.join(output, dataset + "_" + method + "_report.html")
        print("The output path is a directory, adding report file " + dataset + "_" + method + "_report.html to report path.")
    else:
        report_path = output.replace(os.path.splitext(output)[-1], "_" + method + "_report.html")

    if not os.path.exists(os.path.dirname(report_path)):
        os.makedirs(os.path.dirname(report_path))
    
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
    
    size = os.path.getsize(path)
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