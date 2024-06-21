import os
import hashlib
from pathlib import Path
# import tempfile
import shutil
from benchmarks.clustering import clustering_task
from utils.redislogger import *
from utils.mongodb import upsert_benchmarks, upsert_jobs
from utils.unzip import unzip_file_if_compressed
from tools.formating.formating import convert_seurat_sce_to_anndata, load_anndata
from tools.utils.datasplit import sc_train_val_test_split
from fastapi import HTTPException, status
import json
from exceptions.custom_exceptions import CeleryTaskException
from datetime import datetime

def run_data_split(job_id, data_dict:dict):
    datasetId = data_dict['datasetId']
    benchmarksId = data_dict['benchmarksId']
    userID = data_dict['userID']
    adata_path = data_dict['adata_path']
    train_fraction = data_dict['train_fraction']
    validation_fraction = data_dict['validation_fraction']
    test_fraction = data_dict['test_fraction']
    # Serializing the dictionary to a JSON string and encoding to bytes
    encoded_data = json.dumps(data_dict, sort_keys=True).encode('utf-8')
    split_id = hashlib.md5(encoded_data).hexdigest()

    upsert_jobs(
        {
            "job_id": job_id, 
            "created_by": userID,
            "status": "Processing"
        }
    )

    adata_path = unzip_file_if_compressed(job_id, adata_path) 
    # Extract directory and filename from the data filepath
    data_directory = Path(adata_path).parent
    data_filename = Path(adata_path).stem
    adata_dir = Path(f"{data_directory}/{split_id}/")
    train_path = adata_dir / f"{data_filename}_train.h5ad"
    val_path = adata_dir / f"{data_filename}_validation.h5ad"
    test_path = adata_dir / f"{data_filename}_test.h5ad"

    # Ensure the directory exists - added exist_ok=True to prevent any error if directory exists
    os.makedirs(adata_dir, exist_ok=True)

    try:
        adata = load_anndata(adata_path)
        if adata is not None:
            train, validation, test = sc_train_val_test_split(adata, train_fraction, validation_fraction, test_fraction)
        else:
            detail = f'File does not exist at {adata_path}'
            raise CeleryTaskException(detail)
       
        # Define a temporary directory to store the files
        # temp_dir = tempfile.TemporaryDirectory(dir=data_directory)

        # Write AnnData objects to files with unique filenames in the temporary directory
        train.write(str(train_path), compression='gzip')
        validation.write(str(val_path), compression='gzip')
        test.write(str(test_path), compression='gzip')

        # Compress files into a single archive in the same directory
        shutil.make_archive(str(data_directory / f"{data_filename}_data_split"), 'zip', str(adata_dir))

        # Return the path to the compressed archive
        archive_path = str(data_directory / f"{data_filename}_data_split.zip")

        # Updating records using string paths
        upsert_benchmarks(benchmarksId, {
            "datasetId": datasetId,
            "archive_path": archive_path,
            "train_path": str(train_path), 
            "validation_path": str(val_path),
            "test_path": str(test_path)
        })        
        results = {
            "datasetId": datasetId,
            "benchmarksId": benchmarksId,
            "archive_path": archive_path,
            "train_path": str(train_path), 
            "validation_path": str(val_path),
            "test_path": str(test_path)
        }
        
        upsert_jobs(
            {
                "job_id": job_id, 
                "datasetId": datasetId,
                "benchmarksId": benchmarksId,
                "output": archive_path,
                "archive_path": str(archive_path),
                "train_path": str(train_path), 
                "validation_path": str(val_path),
                "test_path": str(test_path),
                "completed_on": datetime.now(),
                "results": results,
                "status": "Success"
            }
        )
        
        return results
    except Exception as e:
        # Handle any errors
        detail=f"Data split is failed: {str(e)}"
        upsert_jobs(
            {
                "job_id": job_id, 
                "results": detail,
                "completed_on": datetime.now(),
                "status": "Failure"
            }
        )
        raise CeleryTaskException(detail)