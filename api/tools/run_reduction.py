import os
from tools.formating.formating import *
from config.celery_utils import get_input_path, get_output
from utils.redislogger import *
from tools.reduction.reduction import run_dimension_reduction, run_clustering
from utils.mongodb import generate_process_id, pp_result_exists, create_pp_results, upsert_jobs
from utils.unzip import unzip_file_if_compressed
from fastapi import HTTPException, status
from exceptions.custom_exceptions import CeleryTaskException
from datetime import datetime
    

def run_reduction(job_id, ds:dict, show_error=True, random_state=0):
    pp_stage = "Summarized"
    process = "Reduction"
    dataset = ds['dataset']
    input = ds['input']
    userID = ds['userID']
    output = ds['output']
    datasetId = ds['datasetId']
    parameters = ds['reduction_params']
    layer = None
    layers = None
    if(len(parameters['layer'].strip())):
        layer = parameters['layer']
    n_neighbors = parameters['n_neighbors']
    n_pcs = parameters['n_pcs']
    resolution = parameters['resolution']
    method = 'UMAP&t-SNE'

    upsert_jobs(
        {
            "job_id": job_id, 
            "created_by": userID,
            "status": "Processing"
        }
    )

    redislogger.info(job_id, f"Using Visualization Parameters: {parameters}")
    
    #Get the absolute path for the given input
    # input = get_input_path(input, userID)
    #Get the absolute path for the given output
    input = unzip_file_if_compressed(job_id, ds['input'])
    md5 = get_md5(input)
    # output = get_output(output, userID, job_id)
    adata = load_anndata(input)
    process_id = generate_process_id(md5, process, method, parameters)
    reduction_results = pp_result_exists(process_id)
    if reduction_results is not None:
        redislogger.info(job_id, "Found existing pre-process results in database, skip dimension reduction.")
    else:
        try:
            redislogger.info(job_id, "Computing PCA, neighborhood graph, tSNE, UMAP, and 3D UMAP")
            adata, msg = run_dimension_reduction(adata, layer=layer, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=random_state)
            if msg is not None: redislogger.warning(job_id, msg)

            redislogger.info(job_id, "Clustering the neighborhood graph.")
            adata = run_clustering(adata, layer=layer, resolution=resolution, random_state=random_state)

            redislogger.info(job_id, "Retrieving metadata and embeddings from AnnData object.")
            reduction_results = get_metadata_from_anndata(adata, pp_stage, process_id, process, method, parameters,  md5, adata_path=output)
            output = get_output_path(output, process_id=process_id, dataset=dataset, method='UMAP_t-SNE')
            adata.write_h5ad(output, compression='gzip')
            adata = None
            reduction_results['datasetId'] = datasetId
            create_pp_results(process_id, reduction_results)  # Insert pre-process results to database
            redislogger.info(job_id, "AnnData object for UMAP & t-SNE reduction is saved successfully")
        except Exception as e:
            # redislogger.error(job_id, "UMAP reduction is failed.")
            detail = f"UMAP or t-SNE reduction is failed: {e}"
            upsert_jobs(
                {
                    "job_id": job_id, 
                    "results": detail,
                    "completed_on": datetime.now(),
                    "status": "Failure"
                }
            )
            os.remove(output)
            raise CeleryTaskException(detail)
    
    if 'layers' in reduction_results.keys():
        layers = reduction_results['layers']

    results = {
        "output": [{method: output}],
        "layers": layers,
        "md5": md5,
        "process_ids": [process_id]
    }

    upsert_jobs(
        {
            "job_id": job_id, 
            "datasetId": datasetId,
            "output": [{method: output}],
            "process_ids": [process_id],
            "layers": layers,
            "results": results,
            "completed_on": datetime.now(),
            "status": "Success"
        }
    )

    return results

        
            

