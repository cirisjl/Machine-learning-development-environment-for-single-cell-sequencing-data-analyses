import pandas as pd
import subprocess
import csv
import os
import scanpy as sc
import gzip

def detect_delimiter(file_path):
    with open(file_path, 'r') as file:
        # Read the first line of the file to detect the delimiter
        first_line = file.readline()
        dialect = csv.Sniffer().sniff(first_line)
        return dialect.delimiter
    
    
def read_text_replace_invalid(file_path, delimiter):
    df = pd.read_csv(file_path, delimiter=delimiter, on_bad_lines='skip', index_col=0)
    df = df.apply(pd.to_numeric, errors='coerce')
    return sc.AnnData(df)


def read_text(file_path):
    if file_path.endswith(".gz"):
        with gzip.open(file_path, 'rt', encoding='utf-8') as file:
            df = pd.read_csv(file, on_bad_lines='skip', index_col=0)
        return sc.AnnData(df)
    else:
        delimiter = detect_delimiter(file_path)
        df = pd.read_csv(file_path, delimiter=delimiter, on_bad_lines='skip', index_col=0)
        return sc.AnnData(df)
    
def load_annData(path, replace_invalid=False):
    show_error=True
    dataset = None
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
            delimiter = detect_delimiter(path)
            if replace_invalid:
                adata = read_text_replace_invalid(path, delimiter)
                print(adata)
                print(adata.var_names[:10])
                print(adata.obs_names[:10])
            else:
                adata = sc.read_text(path, delimiter=detect_delimiter(path))      
        elif path.endswith(".txt.gz"):
            if replace_invalid:
                with gzip.open(path, 'rt') as file:
                    df = pd.read_csv(file, on_bad_lines='skip', index_col=0)
                    df = df.apply(pd.to_numeric, errors='coerce')
                    adata = sc.AnnData(df)
            else:
                adata = sc.read_text(path)      
        elif path.endswith(".gz"):
            adata = sc.read_umi_tools(path)
        elif path.endswith(".h5Seurat") or path.endswith(".h5seurat") or path.endswith(".rds"):
            try:
                current_file = os.path.abspath(__file__)
                # Construct the relative path to the desired file
                relative_path = os.path.join(os.path.dirname(current_file), 'convert_to_anndata.Rmd')

                # Get the absolute path of the desired file
                operation_path = os.path.abspath(relative_path)
                report_path = os.path.join(os.path.dirname(path), "file_conversion_report.html")
                adata_path = os.path.splitext(path)[0] + '.h5ad'
                
                if os.path.exists(adata_path):
                    adata = sc.read_h5ad(adata_path)
                else:
                    s = subprocess.call(["R -e \"rmarkdown::render('" + operation_path + "', params=list(path='" + str(path) + "'), output_file='" + report_path + "')\""], shell = True)
                    print(s)
                    adata = sc.read_h5ad(adata_path)

            except Exception as e:
                print("Object format conversion is failed")
                if show_error: print(e)

    return adata

# def load_annData(file_path, replace_invalid=False):
#     adata = None
#     if os.path.isdir(file_path) and os.path.exists(os.path.join(file_path, "matrix.mtx")) and os.path.exists(
#         os.path.join(file_path, "genes.tsv")) and os.path.exists(os.path.join(file_path, "barcodes.tsv")):
#         adata = sc.read_10x_mtx(file_path,
#                             var_names='gene_symbols',  # use gene symbols for the variable names (variables-axis index)
#                             cache=True)  # write a cache file for faster subsequent reading
#     elif(os.path.exists(file_path)):
#         suffix = os.path.splitext(file_path)[-1]
#         if suffix == ".h5ad":
#             adata = sc.read_h5ad(file_path)
#         elif suffix == ".csv" or suffix == ".tsv":
#             print("Inside the loadAnndata CSV")
#             print(detect_delimiter(file_path))
#             print("Inside the loadAnndata CSV 2")
#             adata = sc.read_csv(file_path, delimiter=detect_delimiter(file_path))
#             print("Inside the loadAnndata CSV 3")
#         elif suffix == ".xlsx" or suffix == ".xls":
#             adata = sc.read_excel(file_path, 0)
#         elif suffix == ".h5" and "pbmc" in file_path:
#             adata = sc.read_10x_h5(file_path)
#         elif suffix == ".h5":
#             adata = sc.read_hdf(file_path)
#         elif suffix == ".loom":
#             adata = sc.read_loom(file_path)
#         elif suffix == ".mtx":
#             adata = sc.read_mtx(file_path)
#         elif suffix == ".txt" or suffix == ".tab" or suffix == ".data":
#             delimiter = detect_delimiter(file_path)
#             if replace_invalid:
#                 adata = read_text_replace_invalid(file_path, delimiter)
#                 print(adata)
#                 print(adata.var_names[:10])
#                 print(adata.obs_names[:10])
#             else:
#                 adata = sc.read_text(file_path, delimiter=detect_delimiter(file_path))
#         elif suffix == ".gz":
#             adata = sc.read_umi_tools(file_path)

#     return adata


def load_invalid_adata(file_path, replace_nan):
    delimiter = detect_delimiter(file_path)
    df = pd.read_csv(file_path, delimiter=delimiter, on_bad_lines='skip', index_col=0)
    if replace_nan == "yes":
        df = df.apply(pd.to_numeric, errors='coerce')
    invalid_rows = df.apply(pd.to_numeric, errors='coerce').isnull().any(axis=1)
    invalid_columns = df.columns[df.apply(pd.to_numeric, errors='coerce').isnull().any()]
    invalid_df = df.loc[invalid_rows, invalid_columns]
    return sc.AnnData(invalid_df)