import os
from tools.formating.formating import *
from config.celery_utils import get_input_path, get_output
from utils.redislogger import *
from utils.mongodb import generate_process_id, pp_result_exists, create_pp_results, upsert_task_results
from utils.unzip import unzip_file_if_compressed
from fastapi import HTTPException, status
from exceptions.custom_exceptions import CeleryTaskException

def run_conversion(task_id, ds:dict, show_error=True):
    outputs = []
    process_ids = []
    pp_stage = "Raw"
    process = "Formatting"
    dataset = ds['dataset']
    input = ds['input']
    userID = ds['userID']
    output = ds['output']
    output_format = ds['output_format']
    
    #Get the absolute path for the given input
    # input = get_input_path(input, userID)
    #Get the absolute path for the given output
    input = unzip_file_if_compressed(task_id, ds['input'])
    # output = get_output(output, userID, task_id)
    adata_path = get_output_path(output, dataset=dataset)
    seurat_path = get_output_path(output, dataset=dataset, format='Seurat')
    sce_path = get_output_path(output, dataset=dataset, format='SingleCellExperiment')
    csv_path = get_output_path(output, dataset=dataset, format='CSV')

    if output_format == "AnnData":
        if os.path.exists(adata_path):
            outputs.append({'AnnData': adata_path})
            redislogger.info(task_id, "AnnData file already exists.")
        else:
            try:
                adata = load_anndata(input) 
                adata.write_h5ad(adata_path, compression='gzip')
                adata = None
                outputs.append({'AnnData': adata_path})
                redislogger.info(task_id, "AnnData object is saved successfully")
            except Exception as e:
                detail = f"Format conversion is failed: {e}"
                os.remove(adata_path)
                raise CeleryTaskException(detail)

    if output_format == "Seurat":
        if os.path.exists(seurat_path):
            outputs.append({'Seurat': seurat_path})
            redislogger.info(task_id, "Seurat file already exists.")
        else:
            try:
                seurat_path = convert_to_seurat_sce(input, seurat_path, format="Seurat") 
                outputs.append({'Seurat': seurat_path})
                redislogger.info(task_id, "Seurat object is saved successfully")
            except Exception as e:
                detail = f"Format conversion is failed: {e}"
                os.remove(seurat_path)
                raise CeleryTaskException(detail)

    if output_format == "SingleCellExperiment":
        if os.path.exists(sce_path):
            outputs.append({'SingleCellExperiment': sce_path})
            redislogger.info(task_id, "SingleCellExperiment file already exists.")
        else:
            try:
                sce_path = convert_to_seurat_sce(input, sce_path, format="SingleCellExperiment") 
                outputs.append({'SingleCellExperiment': sce_path})
                redislogger.info(task_id, "SingleCellExperiment object is saved successfully")
            except Exception as e:
                detail = f"Format conversion is failed: {e}"
                os.remove(sce_path)
                raise CeleryTaskException(detail)

    if output_format == "CSV":
        if os.path.exists(csv_path):
            outputs.append({'CSV': csv_path})
            redislogger.info(task_id, "CSV file already exists.")
        else:
            try:
                # layer = ds['layer']
                adata, counts, csv_path = load_anndata_to_csv(input, csv_path, layer=None)
                adata = None
                outputs.append({'CSV': csv_path})
                redislogger.info(task_id, "CSV file is saved successfully")
            except Exception as e:
                detail = f"Format conversion is failed: {e}"
                os.remove(csv_path)
                raise CeleryTaskException(detail)
        
    results = {
        "taskId": task_id, 
        "owner": userID,
        "inputfile": input,
        "output": outputs,
        "status": "Success"
    }
    upsert_task_results(results)

    return results