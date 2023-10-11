import scanpy as sc
import os
import subprocess
import numpy as np
import pandas as pd
# import anndata2ri
# from rpy2.robjects import r
# anndata2ri.activate()
from detect_delimiter import detect
from string import ascii_letters
import csv


def load_anndata(path, dataset=None, assay='RNA', show_error=True): # assay is optional and only for Seurat object

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
            # print("Inside the loadAnndata CSV 3")
        elif path.endswith(".csv.gz") or path.endswith(".tsv.gz"):
            data = sc.read_csv(path)
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
            adata = sc.read_text(path, delimiter=detect_delim(path))
        elif path.endswith(".txt.gz"):
            adata = sc.read_text(path)
        elif path.endswith(".gz"):
            adata = sc.read_umi_tools(path)
        elif path.endswith(".h5Seurat") or path.endswith(".h5seurat") or path.endswith(".rds"):
            adata_path, assay_names = convert_seurat_sce_to_anndata(path, assay=assay)
            if os.path.exists(adata_path):
                adata = sc.read_h5ad(adata_path)

    return adata


# Convert Seurat/Single-Cell Experiment object to Anndata object and return the path of Anndata object
def convert_seurat_sce_to_anndata(path, assay='RNA'):
    import rpy2.robjects as ro
    # robjects.r.source("formating.R")

    # Get the absolute path of the current file
    current_file = os.path.abspath(__file__)

    # Construct the relative path to the desired file
    relative_path = os.path.join(os.path.dirname(current_file), 'formating.R')

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
            adata_path = str(results[2])
            assay_names = list(results[1])
        except Exception as e:
            print("Object format conversion is failed")
            print(e)

    print("formatting")
    print("AssayNames")
    print(assay_names)
    print("adata_path")
    print(adata_path)
    return adata_path, assay_names


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


def load_anndata_to_csv(input, output, layer, show_error, dataset=None):
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
        output_path = os.path.join(output, dataset + method + ".h5seurat")
        print(output_path)
        print("The output path is a directory, adding output file " + dataset + method + ".h5seurat to the path.")
    
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