import os
from tools.formating.formating import *
from config.celery_utils import get_input_path, get_output
from utils.redislogger import *
from utils.mongodb import generate_process_id, pp_results_exists, create_pp_results, upsert_task_results
from utils.unzip import unzip_file_if_compressed
from fastapi import HTTPException, status
    

def run_conversion(task_id, ds:dict, show_error=True):
    outputs = []
    process_ids = []
    pp_stage = "Raw"
    process = "Formatting"
    dataset = ds['dataset']
    layer = ds['layer']
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
        try:
            adata = load_anndata(input) 
            adata.write_h5ad(adata_path, compression='gzip')
            adata = None
            outputs.append({'adata_path': adata_path})
            redislogger.info(task_id, "AnnData object is saved successfully")
        except Exception as e:
            detail = f"Format conversion is failed: {e}"
            os.remove(adata_path)
            raise HTTPException(
                status_code = status.HTTP_500_INTERNAL_SERVER_ERROR,
                detail = detail
            )

    if output_format == "Seurat":
        try:
            seurat_path = convert_to_seurat_sce(input, seurat_path, format="Seurat") 
            outputs.append({'seurat_path': seurat_path})
            redislogger.info(task_id, "Seurat object is saved successfully")
        except Exception as e:
            detail = f"Format conversion is failed: {e}"
            os.remove(seurat_path)
            raise HTTPException(
                status_code = status.HTTP_500_INTERNAL_SERVER_ERROR,
                detail = detail
            )

    if output_format == "SingleCellExperiment":
        try:
            sce_path = convert_to_seurat_sce(input, sce_path, format="SingleCellExperiment") 
            outputs.append({'sce_path': sce_path})
            redislogger.info(task_id, "SingleCellExperiment object is saved successfully")
        except Exception as e:
            detail = f"Format conversion is failed: {e}"
            os.remove(sce_path)
            raise HTTPException(
                status_code = status.HTTP_500_INTERNAL_SERVER_ERROR,
                detail = detail
            )

    if output_format == "CSV":
        try:
            adata, counts, csv_path = load_anndata_to_csv(input, output, layer=layer)
            adata = None
            csv_path =None
            outputs.append({'csv_path': csv_path})
            redislogger.info(task_id, "CSV file is saved successfully")
        except Exception as e:
            detail = f"Format conversion is failed: {e}"
            os.remove(csv_path)
            raise HTTPException(
                status_code = status.HTTP_500_INTERNAL_SERVER_ERROR,
                detail = detail
            )
        
    results = {
        "taskId": task_id, 
        "owner": userID,
        "inputfile": input,
        "output": outputs,
        "status":"Success"
    }
    upsert_task_results(results)

    return results