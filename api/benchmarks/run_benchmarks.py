from benchmarks.clustering import clustering_task
from utils.redislogger import *
from utils.mongodb import generate_process_id
from utils.unzip import unzip_file_if_compressed
from fastapi import HTTPException, status

def run_benchmarks(task_id, task_dict:dict):
    results = []
    adata_path = task_dict['adata_path']
    label = task_dict['label']
    task_type = task_dict['task_type']
    datasetId = task_dict['datasetId']

    adata_path = unzip_file_if_compressed(adata_path)

    if(task_type=="clustering"):
        try:
            if os.path.exists(adata_path):
                clustering_results = clustering_task(adata_path, label, datasetId, task_id, task_type)
                results.append(clustering_results)
            else:
                raise HTTPException(
                status_code = status.HTTP_500_INTERNAL_SERVER_ERROR,
                detail = f'File does not exist at {adata_path}')
            return results
        except Exception as e:
            # Handle exceptions as needed
            raise HTTPException(status_code=500, detail=f"Clustering benchmarks is failed: {str(e)}")