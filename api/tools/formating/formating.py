import scanpy as sc
import os
import numpy as np
import pandas as pd
# import anndata2ri
# from rpy2.robjects import r
# anndata2ri.activate()
from detect_delimiter import detect
from string import ascii_letters
import csv


def load_anndata(path):

    # path = os.path.abspath(path)
    adata = None
    print(path)

    if os.path.isdir(path) and os.path.exists(os.path.join(path, "matrix.mtx")) and os.path.exists(
            os.path.join(path, "genes.tsv")) and os.path.exists(os.path.join(path, "barcodes.tsv")):
        adata = sc.read_10x_mtx(path,
                             var_names='gene_symbols',  # use gene symbols for the variable names (variables-axis index)
                             cache=True)  # write a cache file for faster subsequent reading
    elif(os.path.exists(path)):
        suffix = os.path.splitext(path)[-1]
        if suffix == ".h5ad":
            adata = sc.read_h5ad(path)
        elif suffix == ".csv" or suffix == ".tsv":
            print("Inside the loadAnndata CSV")
            print(detect_delimiter(path))
            print("Inside the loadAnndata CSV 2")
            adata = sc.read_csv(path, delimiter=detect_delimiter(path))
            print("Inside the loadAnndata CSV 3")
        elif suffix == ".xlsx" or suffix == ".xls":
            adata = sc.read_excel(path, 0)
        elif suffix == ".h5" and "pbmc" in path:
            adata = sc.read_10x_h5(path)
        elif suffix == ".h5":
            adata = sc.read_hdf(path)
        elif suffix == ".loom":
            adata = sc.read_loom(path)
        elif suffix == ".mtx":
            adata = sc.read_mtx(path)
        elif suffix == ".txt" or suffix == ".tab" or suffix == ".data":
            adata = sc.read_text(path, delimiter=detect_delim(path))
        elif suffix == ".gz":
            adata = sc.read_umi_tools(path)

    return adata


def anndata_to_csv(adata, output_path, layer = None):
    counts = None

    if layer is None:
        counts = adata.raw.X.toarray()
    else:
        if type(adata.layers[layer]) != "numpy.ndarray":
            counts = adata.layers[layer].toarray()
        else:
            counts = adata.layers[layer]

    pd.DataFrame(data=counts, index=adata.obs_names, columns=adata.raw.var_names).to_csv(output_path)
    return output_path


def load_anndata_to_csv(input, output, layer, show_error):
    adata = None
    adata_path = None
    counts = None

    if os.path.exists(output):
        try:
            adata = load_anndata(output)
            adata_path = output
        except Exception as e:
            print("File format is not supported.")
            if show_error: print(e)
            return None, None, None
    else:
        try:
            adata = load_anndata(input)
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


def get_output_path(dataset, output, method = '', format = "AnnData"):
    output = os.path.abspath(output)
    if method != '': method = '_' + method
    
    if not os.path.exists(output):
        os.makedirs(output)

    if format == "AnnData":
        output_path = os.path.join(output, dataset + method + ".h5ad")
        print(output_path)
        print("The output path is a directory, adding output file " + dataset + method + ".h5ad to the path.")
    elif format == "SingleCellExperiment":
        output_path = os.path.join(output, dataset + method + ".rds")
        print(output_path)
        print("The output path is a directory, adding output file " + dataset + method + ".rds to the path.")
    elif format == "Seurat":
        output_path = os.path.join(output, dataset + method + ".h5Seurat")
        print(output_path)
        print("The output path is a directory, adding output file " + dataset + method + ".h5Seurat to the path.")
    
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