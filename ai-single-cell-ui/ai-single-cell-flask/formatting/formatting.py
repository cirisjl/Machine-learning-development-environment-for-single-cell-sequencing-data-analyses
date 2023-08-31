import pandas as pd
import csv
import os
import scanpy as sc

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
    delimiter = detect_delimiter(file_path)
    df = pd.read_csv(file_path, delimiter=delimiter, on_bad_lines='skip', index_col=0)
    return sc.AnnData(df)
    

def load_annData(file_path, replace_invalid=False):
    adata = None
    if os.path.isdir(file_path) and os.path.exists(os.path.join(file_path, "matrix.mtx")) and os.path.exists(
        os.path.join(file_path, "genes.tsv")) and os.path.exists(os.path.join(file_path, "barcodes.tsv")):
        adata = sc.read_10x_mtx(file_path,
                            var_names='gene_symbols',  # use gene symbols for the variable names (variables-axis index)
                            cache=True)  # write a cache file for faster subsequent reading
    elif(os.path.exists(file_path)):
        suffix = os.path.splitext(file_path)[-1]
        if suffix == ".h5ad":
            adata = sc.read_h5ad(file_path)
        elif suffix == ".csv" or suffix == ".tsv":
            print("Inside the loadAnndata CSV")
            print(detect_delimiter(file_path))
            print("Inside the loadAnndata CSV 2")
            adata = sc.read_csv(file_path, delimiter=detect_delimiter(file_path))
            print("Inside the loadAnndata CSV 3")
        elif suffix == ".xlsx" or suffix == ".xls":
            adata = sc.read_excel(file_path, 0)
        elif suffix == ".h5" and "pbmc" in file_path:
            adata = sc.read_10x_h5(file_path)
        elif suffix == ".h5":
            adata = sc.read_hdf(file_path)
        elif suffix == ".loom":
            adata = sc.read_loom(file_path)
        elif suffix == ".mtx":
            adata = sc.read_mtx(file_path)
        elif suffix == ".txt" or suffix == ".tab" or suffix == ".data":
            delimiter = detect_delimiter(file_path)
            if replace_invalid:
                adata = read_text_replace_invalid(file_path, delimiter)
                print(adata)
                print(adata.var_names[:10])
                print(adata.obs_names[:10])
            else:
                adata = sc.read_text(file_path, delimiter=detect_delimiter(file_path))
        elif suffix == ".gz":
            adata = sc.read_umi_tools(file_path)

    return adata


def load_invalid_adata(file_path, replace_nan):
    delimiter = detect_delimiter(file_path)
    df = pd.read_csv(file_path, delimiter=delimiter, on_bad_lines='skip', index_col=0)
    if replace_nan == "yes":
        df = df.apply(pd.to_numeric, errors='coerce')
    invalid_rows = df.apply(pd.to_numeric, errors='coerce').isnull().any(axis=1)
    invalid_columns = df.columns[df.apply(pd.to_numeric, errors='coerce').isnull().any()]
    invalid_df = df.loc[invalid_rows, invalid_columns]
    return sc.AnnData(invalid_df)