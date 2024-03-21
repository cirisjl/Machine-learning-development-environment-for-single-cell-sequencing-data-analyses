import os
import subprocess
import sys

from tools.qc.scanpy_qc import run_scanpy_qc
from tools.qc.dropkick_qc import run_dropkick_qc
from tools.qc.scrublet_calls import predict_scrublet
# sys.path.append('..')
from tools.formating.formating import *
from config.celery_utils import get_input_path, get_output


def run_qc(task_id, dataset, input,userID, output, methods, idtype='SYMBOL', colour_by='NULL', shape_by_1='NULL', shape_by_2='NULL', default_assay='RNA', show_error=True):
    if methods is None:
        print("No quality control method is selected.")
        return None   
    

     #Get the absolute path for the given input
    # input = get_input_path(input, userID)
    #Get the absolute path for the given output
    output = get_output(output, userID,task_id)
    methods = [x.upper() for x in methods if isinstance(x,str)]


    
    if "SCANPY" in methods or "DROPKICK" in methods:
        adata = load_anndata(input)

        if adata is None:
            print("File format is not supported.")
            return None

        # Scanpy QC
        if "SCANPY" in methods:
            try:
                adata = run_scanpy_qc(adata, userID)
                output_path = get_output_path(dataset, output, method='scanpy')
                 # Save AnnData object
                adata.write_h5ad(output_path, compression='gzip')
                print("AnnData object for Scanpy QC is saved successfully")
            except Exception as e:
                print("Scanpy QC is failed")
                if show_error: print(e)

        # Dropkick QC
        if "DROPKICK" in methods:
            try:
                adata = run_dropkick_qc(adata)
                output_path = get_output_path(dataset, output, method='dropkick')
                 # Save AnnData object
                adata.write_h5ad(output_path, compression='gzip')
                print("AnnData object for Dropkick QC is saved successfully")
            except Exception as e:
                print("Dropkick QC is failed")
                if show_error: print(e)

    # Bioconductor QC
    if "BIOCONDUCTOR" in methods:
        try:
            output_path = get_output_path(dataset, output, method='Bioconductor', format='SingleCellExperiment')
            report_path = get_report_path(dataset, output_path, "Bioconductor")

            # Get the absolute path of the current file
            current_file = os.path.abspath(__file__)

            # Construct the relative path to the desired file
            relative_path = os.path.join(os.path.dirname(current_file), 'qc', 'bioconductor_qc.Rmd')

            # Get the absolute path of the desired file
            bioconductor_path = os.path.abspath(relative_path)
            
            # bioconductor_path = os.path.abspath("qc/bioconductor_qc.Rmd")
            s = subprocess.call(["R -e \"rmarkdown::render('" + bioconductor_path + "', params=list(dataset='" + str(dataset) + "', input_path='" + input + "', idtype='" + idtype + "', colour_by='" + colour_by + "', shape_by_1='" + shape_by_1 + "', shape_by_2='" + shape_by_2 + "', output='" + output_path + "', output_format='SingleCellExperiment'), output_file='" + report_path + "')\""], shell = True)
            print(s)
        except Exception as e:
            print("Bioconductor QC is failed")
            if show_error: print(e)

    # Seurat QC
    if "SEURAT" in methods:
        try:
            path_of_scrublet_calls = predict_scrublet(input)
            output_path = get_output_path(dataset, output, method='Seurat', format='Seurat')
            report_path = get_report_path(dataset, output_path, "Seurat")
            
            # Get the absolute path of the current file
            current_file = os.path.abspath(__file__)

            # Construct the relative path to the desired file
            relative_path = os.path.join(os.path.dirname(current_file), 'qc', 'seurat_qc.Rmd')

            # Get the absolute path of the desired file
            seurat_path = os.path.abspath(relative_path)
            # seurat_path = os.path.abspath("seurat_qc.Rmd")
            s = subprocess.call(["R -e \"rmarkdown::render('" + seurat_path + "', params=list(dataset='" + str(dataset) + "', input='" + input + "', default_assay='" + default_assay + "', output='" + output_path + "', output_format='Seurat', path_of_scrublet_calls='" + path_of_scrublet_calls + "'), output_file='" + report_path + "')\""], shell = True)
            print(s)
        except Exception as e:
            print("Seurat QC is failed")
            if show_error: print(e)

    return {'status': 'Success'}