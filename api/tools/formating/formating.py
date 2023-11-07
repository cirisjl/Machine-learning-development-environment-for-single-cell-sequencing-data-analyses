import scanpy as sc
import os
import sys
import subprocess
import numpy as np
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
# anndata2ri.activate()
from detect_delimiter import detect
from string import ascii_letters
import csv
import gzip

# Defining the R script and loading the instance in Python
ro.r['source'](os.path.abspath(os.path.join(os.path.dirname(__file__), 'formating.R')))


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
            adata_path, assay_names = convert_seurat_sce_to_anndata(path, assay=assay)
            if os.path.exists(adata_path):
                adata = sc.read_h5ad(adata_path)

    return adata


def change_file_extension(file_path, new_extension):
    # Split the file path into directory and filename
    directory, base_filename = os.path.split(file_path)

    # Split the base filename into name and current extension
    name, current_extension = os.path.splitext(base_filename)

    # Create the new file path with the desired extension
    new_file_path = os.path.join(directory, f"{name}.{new_extension}")

    return new_file_path

def get_metadata_from_seurat(path):
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
    except Exception as e:
        print(e)

    return default_assay, assay_names, metadata, nCells, nGenes, genes, cells, HVGsID, pca, tsne, umap


# Convert Seurat/Single-Cell Experiment object to Anndata object and return the path of Anndata object
def convert_seurat_sce_to_anndata(path, assay='RNA'):
    # Access the loaded R functions
    ConvertSeuratSCEtoAnndata_r = ro.globalenv['ConvertSeuratSCEtoAnndata']

    assay_names = None
    adata_path = None

    if path.endswith(".h5Seurat") or path.endswith(".h5seurat") or path.endswith(".rds") or path.endswith(".Robj"):
        try:
            results = ConvertSeuratSCEtoAnndata_r(path, assay=assay)
            adata_path = str(results.rx2('anndata_path'))
            assay_names = list(results[1])
        except Exception as e:
            print("Object format conversion is failed")
            print(e)

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