import os
from tools.formating.formating import *
from config.celery_utils import get_input_path, get_output
from utils.redislogger import *
from tools.reduction.reduction import run_dimension_reduction, run_clustering
from utils.mongodb import generate_process_id, pp_results_exists, create_pp_results
from utils.unzip import unzip_file_if_compressed
    

def run_reduction(task_id, ds:dict, show_error=True, random_state=0):
    results = []
    pp_stage = "Summarized"
    process = "Reduction"
    dataset = ds['dataset']
    layer = ds['layer']
    input = ds['input']
    userID = ds['userID']
    output = ds['output']
    methods = "UMAP"
    parameters = ds['reduction_params']
    n_neighbors = parameters['n_neighbors']
    n_pcs = parameters['n_pcs']
    resolution = parameters['resolution']
    
    #Get the absolute path for the given input
    # input = get_input_path(input, userID)
    #Get the absolute path for the given output
    input = unzip_file_if_compressed(input)
    md5 = get_md5(input)
    output = get_output(output, userID, task_id)

    adata = load_anndata(input)
    method='MAGIC'
    process_id = generate_process_id(md5, process, method, parameters)
    reduction_results = pp_results_exists(process_id)
    if reduction_results is not None:
        redislogger.info(task_id, "Found existing pre-process results in database, skip dimension reduction.")
    else:
        try:
            redislogger.info(task_id, "Computing PCA, neighborhood graph, tSNE, UMAP, and 3D UMAP")
            adata, msg = run_dimension_reduction(adata, layer=layer, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=random_state)
            if msg is not None: redislogger.warning(task_id, msg)

            redislogger.info(task_id, "Clustering the neighborhood graph.")
            adata = run_clustering(adata, layer=layer, resolution=resolution, random_state=random_state)

            redislogger.info(task_id, "Retrieving metadata and embeddings from AnnData object.")
            reduction_results = get_metadata_from_anndata(adata, pp_stage, process_id, process, method, parameters, adata_path=output)

            output = get_output_path(dataset, output, method='UMAP')
            adata.write_h5ad(output, compression='gzip')
            adata = None
            create_pp_results(reduction_results)  # Insert pre-process results to database
            redislogger.info(task_id, "AnnData object for UMAP reduction is saved successfully")
        except Exception as e:
            redislogger.error(task_id, "UMAP reduction is failed.")
            if show_error: redislogger.error(task_id, f"UMAP reduction is failed: {e}")
        
    results.append({
        "task_id": task_id, 
        "inputfile": input,
        "layers": reduction_results.layers,
        "md5": md5,
        "process_id": process_id,
        "pp_results": reduction_results,
        "message": "Dimension reduction completed successfully."
    })

    return results

        
            

