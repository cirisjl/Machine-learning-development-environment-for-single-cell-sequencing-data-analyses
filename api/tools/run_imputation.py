import os
import subprocess
import sys
from tools.formating.formating import *
from tools.imputation.MAGIC import magic_impute
from config.celery_utils import get_input_path, get_output
from utils.redislogger import *
from tools.reduction.reduction import run_dimension_reduction, run_clustering
from utils.mongodb import generate_process_id, pp_results_exists, create_pp_results
from utils.unzip import unzip_file_if_compressed
from fastapi import HTTPException, status
    

def run_imputation(task_id, ds:dict, show_error=True, random_state=0):
    results = []
    pp_results = []
    process_ids = []
    pp_stage = "Corrected"
    process = "Imputation"
    dataset = ds['dataset']
    layer = ds['layer']
    input = ds['input']
    userID = ds['userID']
    output = ds['output']
    methods = ds['methods']
    parameters = ds['imputation_params']
    genes = parameters['genes']
    ncores = parameters['ncores']
    n_neighbors = parameters['n_neighbors']
    n_pcs = parameters['n_pcs']
    resolution = parameters['resolution']
    status = 'Successful'
    
    if methods is None:
        redislogger.warning(task_id, "No imputation method is selected.")
        return None
    
    #Get the absolute path for the given input
    # input = get_input_path(input, userID)
    #Get the absolute path for the given output
    input = unzip_file_if_compressed(task_id, ds['input'])
    md5 = get_md5(input)
    output = get_output(output, userID, task_id)
    methods = [x.upper() for x in methods if isinstance(x,str)]
    
    if "MAGIC" in methods:
        adata = load_anndata(input)
        method='MAGIC'
        process_id = generate_process_id(md5, process, method, parameters)
        imputation_results = pp_results_exists(process_id)
        if imputation_results is not None:
            redislogger.info(task_id, "Found existing pre-process results in database, skip MAGIC imputation.")
        else:
            if adata is None:
                redislogger.error(task_id, f"File format is not supported: {input}")
            elif 'MAGIC' not in adata.layers.keys(): 
                try:
                    redislogger.info(task_id, "Start Magic imputation...")
                    counts = adata.X
                    data_magic = magic_impute(counts, genes)
                    adata.layers['MAGIC'] = data_magic

                    redislogger.info(task_id, "Computing PCA, neighborhood graph, tSNE, UMAP, and 3D UMAP")
                    adata, msg = run_dimension_reduction(adata, layer='MAGIC', n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=random_state)
                    if msg is not None: redislogger.warning(task_id, msg)

                    redislogger.info(task_id, "Clustering the neighborhood graph.")
                    adata = run_clustering(adata, layer='MAGIC', resolution=resolution, random_state=random_state)

                    redislogger.info(task_id, "Retrieving metadata and embeddings from AnnData object.")
                    imputation_results = get_metadata_from_anndata(adata, pp_stage, process_id, process, method, parameters, adata_path=output)

                    output = get_output_path(dataset, output, method='MAGIC_imputation')
                    adata.write_h5ad(output, compression='gzip')
                    adata = None
                    redislogger.info(task_id, "AnnData object for MAGIC imputation is saved successfully")
                except Exception as e:
                    detail = f"MAGIC imputation is failed: {e}"
                    raise HTTPException(
                        status_code = status.HTTP_500_INTERNAL_SERVER_ERROR,
                        detail = detail
                    )
            else: 
                redislogger.warning(task_id, "'MAGIC_imputed' layer already exists.")
                imputation_results = get_metadata_from_anndata(adata, pp_stage, process_id, process, method, parameters, adata_path=output)

            create_pp_results(imputation_results)  # Insert pre-process results to database

        pp_results.append(imputation_results)
        process_ids.append(imputation_results)
        

    # if "scGNN" in methods:
    #     if 'scGNN_imputed' not in adata.layers.keys(): 
    #         try:
    #             output = get_output_path(dataset, output, method='scGNN_imputation')
    #             redislogger.info(task_id, "AnnData object for scGNN imputation is saved successfully")          
    #         except Exception as e:
    #             redislogger.error(task_id, "scGNN imputation is failed.")
    #             if show_error: redislogger.error(task_id, f"scGNN imputation is failed: {e}")
    #     else: 
    #         redislogger.warning(task_id, "'scGNN_imputed' layer already exists.") 
    

    if "SAVER" in methods:
        method='SAVER'
        process_id = generate_process_id(md5, process, method, parameters)
        imputation_results = pp_results_exists(process_id)

        if imputation_results is not None:
            redislogger.info(task_id, "Found existing pre-process results in database, skip SAVER imputation.")
        else:
            adata, counts, csv_path = load_anndata_to_csv(input, output, layer, show_error)
            if adata is None:
                redislogger.warning(task_id, f"File format is not supported: {input}")
            elif 'SAVER' not in adata.layers.keys(): 
                try:
                    output = get_output_path(dataset, output, method='SAVER_imputation')
                    report_path = get_report_path(dataset, output, "SAVER")
                    
                    # Get the absolute path of the current file
                    current_file = os.path.abspath(__file__)

                    # Construct the relative path to the desired file
                    relative_path = os.path.join(os.path.dirname(current_file), 'imputation', 'SAVER.Rmd')

                    # Get the absolute path of the desired file
                    saver_path = os.path.abspath(relative_path)

                    # saver_path = os.path.abspath("imputation/SAVER.Rmd")
                    s = subprocess.call(["R -e \"rmarkdown::render('" + saver_path + "', params=list(dataset='" + str(dataset) + "', input='" + csv_path + "', output='" + output + "', output_format='AnnData', ncores=" + str(ncores) + "), output_file='" + report_path + "')\""], shell = True)
                    redislogger.info(task_id, s)

                    if os.path.exists(output):
                        adata = load_anndata(output)
                        redislogger.info(task_id, "Computing PCA, neighborhood graph, tSNE, UMAP, and 3D UMAP")
                        adata, msg = run_dimension_reduction(adata, layer='SAVER', n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=random_state)
                        if msg is not None: redislogger.warning(task_id, msg)

                        redislogger.info(task_id, "Clustering the neighborhood graph.")
                        adata = run_clustering(adata, layer='SAVER', resolution=resolution, random_state=random_state)

                        redislogger.info(task_id, "Retrieving metadata and embeddings from AnnData object.")
                        imputation_results = get_metadata_from_anndata(adata, pp_stage, process_id, process, method, parameters, adata_path=output)

                        adata.write_h5ad(output, compression='gzip')
                        adata = None
                        redislogger.info(task_id, "AnnData object for SAVER imputation is saved successfully")
                    else:
                        raise ValueError("AnnData file does not exist due to the failure of Bioconductor QC.")
            
                except Exception as e:
                    detail = f"SAVER imputation is failed: {e}"
                    raise HTTPException(
                        status_code = status.HTTP_500_INTERNAL_SERVER_ERROR,
                        detail = detail
                    )
            else: 
                redislogger.warning(task_id, "'SAVER' layer already exists.")
                imputation_results = get_metadata_from_anndata(adata, pp_stage, process_id, process, method, parameters, adata_path=output)
            create_pp_results(imputation_results)  # Insert pre-process results to database

        pp_results.append(imputation_results)
        process_ids.append(imputation_results)
        
    results.append({
        "task_id": task_id, 
        "inputfile": input,
        "md5": md5,
        "process_id": process_ids,
        "pp_results": pp_results,
        "status":"Imputation completed successfully."
    })

    return results

    
    

        
            

