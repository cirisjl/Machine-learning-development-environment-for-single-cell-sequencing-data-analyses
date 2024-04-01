import os
import hashlib
from pathlib import Path
import tempfile
import shutil
from benchmarks.clustering import clustering_task
from utils.redislogger import *
from utils.mongodb import generate_process_id
from utils.unzip import unzip_file_if_compressed
from tools.formating.formating import convert_seurat_sce_to_anndata, load_anndata
from tools.utils.datasplit import subset_by_obskey
from fastapi import HTTPException, status


def run_subset_data(task_id, data_dict:dict):
    adata_path = data_dict['adata_path']
    obskey = data_dict['obskey']
    values = data_dict['values']
    subset_id = hashlib.md5(f"{data_dict}".encode("utf_8")).hexdigest()
    adata_path = unzip_file_if_compressed(adata_path)
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
            raise HTTPException(
            status_code = status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail = f'File does not exist at {adata_path}')

       
        adata_sub.write(archive_path, compression='gzip')

        return {"archive_path": archive_path,
                "status": "AnnData is subset successfully.",}
    
    except Exception as e:
        # Handle any errors
        raise HTTPException(status_code=500, detail=f"Subsetting data is failed: {str(e)}")
