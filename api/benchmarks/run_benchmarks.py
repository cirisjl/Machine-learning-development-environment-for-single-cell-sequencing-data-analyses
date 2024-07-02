from benchmarks.clustering import clustering_task
from utils.redislogger import *
from utils.mongodb import upsert_benchmarks, upsert_jobs
from utils.unzip import unzip_file_if_compressed
from fastapi import HTTPException, status
from exceptions.custom_exceptions import CeleryTaskException
from datetime import datetime

def run_benchmarks(job_id, task_dict:dict):
    adata_path = task_dict['adata_path']
    label = task_dict['label']
    task_type = task_dict['task_type']
    benchmarksId = task_dict['benchmarksId']
    datasetId = task_dict['datasetId']
    userID = task_dict['userID']

    upsert_jobs(
        {
            "job_id": job_id, 
            "created_by": userID,
            "status": "Processing"
        }
    )

    adata_path = unzip_file_if_compressed(job_id,adata_path)
    
    if(task_type=="Clustering"):
        try:
            if os.path.exists(adata_path):
                clustering_results = clustering_task(adata_path, label, benchmarksId, datasetId, job_id, task_type)
                upsert_benchmarks(benchmarksId, clustering_results)
                results = {
                    "datasetId": datasetId,
                    "benchmarksId": benchmarksId
                }

                upsert_jobs(
                    {
                        "job_id": job_id, 
                        "datasetId": datasetId,
                        "benchmarksId": benchmarksId,
                        "results": clustering_results,
                        "completed_on": datetime.now(),
                        "status": "Success"
                    }
                )

                return clustering_results
            else:
                detail = f'File does not exist at {adata_path}'
                raise CeleryTaskException(detail)
            
        except Exception as e:
            # Handle exceptions as needed
            detail=f"Clustering benchmarks is failed: {str(e)}"
            upsert_jobs(
                {
                    "job_id": job_id, 
                    "results": detail,
                    "completed_on": datetime.now(),
                    "status": "Failure"
                }
            )
            raise CeleryTaskException(detail)