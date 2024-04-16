from benchmarks.clustering import clustering_task
from utils.redislogger import *
from utils.mongodb import upsert_benchmarks, upsert_task_results
from utils.unzip import unzip_file_if_compressed
from fastapi import HTTPException, status

def run_benchmarks(task_id, task_dict:dict):
    adata_path = task_dict['adata_path']
    label = task_dict['label']
    task_type = task_dict['task_type']
    benchmarksId = task_dict['benchmarksId']
    datasetId = task_dict['datasetId']
    userID = task_dict['userID']

    adata_path = unzip_file_if_compressed(task_id,adata_path)

    if(task_type=="Clustering"):
        try:
            if os.path.exists(adata_path):
                clustering_results = clustering_task(adata_path, label, datasetId, task_id, task_type)
                upsert_benchmarks(benchmarksId, clustering_results)
                results = {
                    "taskId": task_id,
                    "owner": userID,
                    "datasetId": datasetId,
                    "benchmarksId": benchmarksId,
                    "status": "Success"
                }
                upsert_task_results(results)
                return results
            else:
                raise HTTPException(
                status_code = status.HTTP_500_INTERNAL_SERVER_ERROR,
                detail = f'File does not exist at {adata_path}')
            
        except Exception as e:
            # Handle exceptions as needed
            raise HTTPException(status_code=500, detail=f"Clustering benchmarks is failed: {str(e)}")