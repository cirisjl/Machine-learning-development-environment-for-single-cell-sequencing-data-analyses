import pandas as pd
import subprocess
import csv
import os
import scanpy as sc
from detect_delimiter import detect
import gzip

def detect_delimiter(file_path):
    with open(file_path, 'r') as file:
        # Read the first line of the file to detect the delimiter
        first_line = file.readline()
        dialect = csv.Sniffer().sniff(first_line)
        return dialect.delimiter
    
def convert_gz_to_txt(gz_file_path, txt_file_path):
  with gzip.open(gz_file_path, 'rb') as f_in:
    with open(txt_file_path, 'w') as f_out:
      f_out.write(f_in.read().decode())
    
def detect_delim(path):
    # look at the first ten thousand bytes to guess the character encoding
    with open(path, 'rb') as file:
        rawdata = file.read(10000)
        rawdata = rawdata.decode('utf-8')
        delimiter = detect(rawdata, whitelist=[' ', ',', ';', ':', '|', '\t'])
        return delimiter
    
    
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
    # file_path = '/usr/src/app/storage/kbcfh/GZ Dataset/GSE50244_Genes_counts_TMM_NormLength_atLeastMAF5_expressed.txt.gz'
    if file_path.endswith(".gz"):
        # with gzip.open(file_path, 'rb') as file:
        #     compressed_data = file.read()
       
        # decompressed_data = gzip.decompress(compressed_data)
        # # Get the base name of the file without extension
        
        file_name_without_extension = os.path.splitext(os.path.basename(file_path))[0]

        # Create the new file path with a .txt extension
        new_file_path = os.path.join(os.path.dirname(file_path), file_name_without_extension + '.txt')
        convert_gz_to_txt(file_path, new_file_path)

        # Write the decompressed data to a plain .txt file
        # with open(new_file_path, 'wb') as txt_file:
        #     txt_file.write(decompressed_data)

        df = pd.read_csv(new_file_path, sep="\t", on_bad_lines='skip', index_col=0)
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

    #Check for the corrected / updated file by the user. If the file exits then show the contents of the updated file.
    file_name = path.split("/")
    fileparts = file_name[len(file_name)-1].split(".")
    filename = fileparts[0] + "_user_corrected.h5ad"
    updated_filename = os.path.join(os.path.dirname(path), filename)
    if os.path.exists(updated_filename):
        path = updated_filename

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
                adata = read_text_replace_invalid(path, "/t")
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


# Convert Seurat/Single-Cell Experiment object to Anndata object and return the path of Anndata object
def convert_seurat_sce_to_anndata(path, assay='RNA'):
    import rpy2.robjects as ro

    # Get the absolute path of the current file
    current_file = os.path.abspath(__file__)

    # Construct the relative path to the desired file
    relative_path = os.path.join(os.path.dirname(current_file), 'formatting.R')

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
            adata_path = str(results.rx2('anndata_path'))
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