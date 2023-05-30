import os
import subprocess
import sys

from tools.qc.scanpy_qc import scanpy_qc
from tools.qc.dropkick_qc import dropkick_qc
# sys.path.append('..')
from tools.formating.formating import *

def qc(dataset, input, output, methods, path_of_scrublet_calls='./scrublet_calls.tsv', show_error=True):
    if methods is None:
        print("No imputation method is selected.")
        return None   
    
    methods = [x.upper() for x in methods if isinstance(x,str)]
    
    if "SCANPY" in methods or "DROPKICK" in methods:
        adata = load_anndata(input)

        if adata is None:
            print("File format is not supported.")
            return None

        # Scanpy QC
        if "SCANPY" in methods:
            try:
                adata = scanpy_qc(adata)
                output = output_path_check(dataset, output, method='scanpy')
                 # Save AnnData object
                adata.write_h5ad(output, compression='gzip')
                print("AnnData object for Scanpy QC is saved successfully")
            except Exception as e:
                print("Scanpy QC is failed")
                if show_error: print(e)

        # Dropkick QC
        if "DROPKICK" in methods:
            try:
                adata = dropkick_qc(adata)
                output = output_path_check(dataset, output, method='dropkick')
                 # Save AnnData object
                adata.write_h5ad(output, compression='gzip')
                print("AnnData object for Dropkick QC is saved successfully")
            except Exception as e:
                print("Dropkick QC is failed")
                if show_error: print(e)

    # Bioconductor QC
    if "BIOCONDUCTOR" in methods:
        try:
            output = output_path_check(dataset, output, method='Bioconductor', format='SingleCellExperiment')
            report_path = get_report_path(dataset, output, "Bioconductor")
            bioconductor_path = os.path.abspath("bioconductor_qc.Rmd")
            s = subprocess.call(["R -e \"rmarkdown::render('" + bioconductor_path + "', params=list(dataset='" + str(dataset) + "', input='" + input + "', output='" + output + "', output_format='SingleCellExperiment'), output_file='" + report_path + "')\""], shell = True)
            print(s)
        except Exception as e:
            print("Bioconductor QC is failed")
            if show_error: print(e)

    # Seurat QC
    if "SEURAT" in methods:
        try:
            output = output_path_check(dataset, output, method='Seurat')
            report_path = get_report_path(dataset, output, "Seurat")
            seurat_path = os.path.abspath("seurat_qc.Rmd")
            s = subprocess.call(["R -e \"rmarkdown::render('" + seurat_path + "', params=list(dataset='" + str(dataset) + "', input='" + input + "', output='" + output + "', output_format='SingleCellExperiment', path_of_scrublet_calls='" + path_of_scrublet_calls + "'), output_file='" + report_path + "')\""], shell = True)
            print(s)
        except Exception as e:
            print("Seurat QC is failed")
            if show_error: print(e)

    return {'status': 'Success'}