import os
import hashlib
from pathlib import Path
import tempfile
import shutil
from benchmarks.clustering import clustering_task
from utils.redislogger import *
from utils.mongodb import upsert_benchmarks, upsert_task_results
from utils.unzip import unzip_file_if_compressed
from tools.formating.formating import load_anndata
from tools.utils.datasplit import subset_by_obskey
from fastapi import HTTPException, status
from exceptions.custom_exceptions import CeleryTaskException


def run_subset_data(task_id, data_dict:dict):
    results = []
    datasetId = data_dict['datasetId']
    benchmarksId = data_dict['benchmarksId']
    userID = data_dict['userID']
    adata_path = data_dict['adata_path']
    obskey = data_dict['obskey']
    values = data_dict['values']
    subset_id = hashlib.md5(f"{data_dict}".encode("utf_8")).hexdigest()
    adata_path = unzip_file_if_compressed(task_id, adata_path)
    adata_sub = None
    # Extract directory and filename from the data filepath
    data_directory = Path(adata_path).parent
    data_filename = Path(adata_path).stem
    archive_dir = f"{data_directory}/{subset_id}/"

    if not os.path.exists(archive_dir):
        os.makedirs(archive_dir)
    
    # Return the path to the compressed archive
    archive_path = f"{archive_dir}/{data_filename}_sub.h5ad"

    try:
        adata = load_anndata(adata_path)
        if adata is not None:
            adata_sub = subset_by_obskey(adata, obskey, values)
        else:
            detail = f'File does not exist at {adata_path}'
            raise CeleryTaskException(detail)
        
        adata_sub.write(archive_path, compression='gzip')
        upsert_benchmarks(benchmarksId, {"adata_path": archive_path})

        results.append(
            {
                "taskId": task_id,
                "owner": userID,
                "datasetId": datasetId,
                "benchmarksId": benchmarksId,
                "adata_path": archive_path,
                "status": "Success"
            }
        )

        upsert_task_results(results)
        return results
    
    except Exception as e:
        # Handle any errors
        detail=f"Subsetting data is failed: {str(e)}"
        raise CeleryTaskException(detail)
