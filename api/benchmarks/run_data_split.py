import os
import hashlib
from pathlib import Path
# import tempfile
import shutil
from benchmarks.clustering import clustering_task
from utils.redislogger import *
from utils.mongodb import generate_process_id
from utils.unzip import unzip_file_if_compressed
from tools.formating.formating import convert_seurat_sce_to_anndata, load_anndata
from tools.utils.datasplit import sc_train_val_test_split
from fastapi import HTTPException, status

def run_data_split(task_id, data_dict:dict):
    adata_path = data_dict['adata_path']
    train_fraction = data_dict['train_fraction']
    validation_fraction = data_dict['validation_fraction']
    test_fraction = data_dict['test_fraction']
    split_id = hashlib.md5(f"{data_dict}").hexdigest()

    adata_path = unzip_file_if_compressed(adata_path) 
    # Extract directory and filename from the data filepath
    data_directory = Path(adata_path).parent
    data_filename = Path(adata_path).stem
    adata_dir = f"{data_directory}/{split_id}/"
    train_path = f"{adata_dir}/{data_filename}_train.h5ad"
    val_path = f"{adata_dir}/{data_filename}_validation.h5ad"
    test_path = f"{adata_dir}/{data_filename}_test.h5ad"
    if not os.path.exists(adata_dir):
        os.makedirs(adata_dir)

    try:
        adata = load_anndata(adata_path)
        if adata is not None:
            train, validation, test = sc_train_val_test_split(adata, train_fraction, validation_fraction, test_fraction)
        else:
            raise HTTPException(
            status_code = status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail = f'File does not exist at {adata_path}')
       
        # Define a temporary directory to store the files
        # temp_dir = tempfile.TemporaryDirectory(dir=data_directory)

        # Write AnnData objects to files with unique filenames in the temporary directory
        train.write(train_path, compression='gzip')
        validation.write(val_path, compression='gzip')
        test.write(test_path, compression='gzip')

        # Compress files into a single archive in the same directory
        shutil.make_archive(data_directory / f"{data_filename}_data_split", 'zip', adata_dir)

        # Return the path to the compressed archive
        archive_path = data_directory / f"{data_filename}_data_split.zip"
        return {"archive_path": archive_path,
                "train_path": train_path, 
                "validation_path": val_path,
                "test_path": test_path, 
                "status": "Data split successfully."}
    except Exception as e:
        # Handle any errors
        raise HTTPException(status_code=500, detail=f"Data split is failed: {str(e)}")