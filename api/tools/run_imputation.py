import os
import subprocess
import sys
from tools.formating.formating import *
from tools.imputation.MAGIC import magic_impute
from config.celery_utils import get_input_path, get_output
from utils.redislogger import *
from tools.reduction.reduction import run_dimension_reduction, run_clustering
from utils.mongodb import generate_process_id, pp_result_exists, create_pp_results, upsert_jobs
from utils.unzip import unzip_file_if_compressed
from fastapi import HTTPException, status
from exceptions.custom_exceptions import CeleryTaskException
from datetime import datetime
    

def run_imputation(job_id, ds:dict, show_error=True, random_state=0):
    pp_results = []
    process_ids = []
    imputation_output = []
    pp_stage = "Corrected"
    process = "Imputation"
    dataset = ds['dataset']
    input = ds['input']
    userID = ds['userID']
    output = ds['output']
    datasetId = ds['dataset_id']
    parameters = ds['imputation_params']
    layer = None
    if parameters['layer'] is not None and parameters['layer'].strip != "":
        layer = parameters['layer']

    methods = parameters['methods']
    genes = parameters['genes']
    ncores = parameters['ncores']
    n_neighbors = parameters['n_neighbors']
    n_pcs = parameters['n_pcs']
    resolution = parameters['resolution']

    upsert_jobs(
        {
            "job_id": job_id, 
            "created_by": userID,
            "status": "Processing"
        }
    )
    
    if methods is None:
        redislogger.error(job_id, "No imputation method is selected.")
        detail = 'No imputation method is selected.'
        raise CeleryTaskException(detail)
    
    #Get the absolute path for the given input
    # input = get_input_path(input, userID)
    #Get the absolute path for the given output
    input = unzip_file_if_compressed(job_id, ds['input'])
    md5 = get_md5(input)
    # output = get_output(output, userID, job_id)

    methods = [x.upper() for x in methods if isinstance(x,str)]
    
    if "MAGIC" in methods:
        adata = None
        method='MAGIC'
        process_id = generate_process_id(md5, process, method, parameters)
        imputation_results = pp_result_exists(process_id)
        output = get_output_path(output, process_id, dataset, method='MAGIC_imputation')
        
        if imputation_results is not None:
            redislogger.info(job_id, "Found existing pre-process results in database, skip MAGIC imputation.")
            process_ids.append(process_id)
            imputation_output = imputation_results["outputs"]

        else:
            if os.path.exists(output): # If output exist from the last run, then just pick up it.
                redislogger.info(job_id, "Output already exists, start from the result of the last run.")
                adata = load_anndata(output)
                redislogger.info(job_id, "Computing PCA, neighborhood graph, tSNE, UMAP, and 3D UMAP")
                adata, msg = run_dimension_reduction(adata, layer='MAGIC', n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=random_state)
                if msg is not None: redislogger.warning(job_id, msg)

                redislogger.info(job_id, "Clustering the neighborhood graph.")
                adata = run_clustering(adata, layer='MAGIC', resolution=resolution, random_state=random_state)

                redislogger.info(job_id, "Retrieving metadata and embeddings from AnnData object.")
                imputation_results = get_metadata_from_anndata(adata, pp_stage, process_id, process, method, parameters, md5, adata_path=output, layer='MAGIC')
                adata.write_h5ad(output, compression='gzip')
                imputation_output.append({"MAGIC": output})
                imputation_results["outputs"] = imputation_output
                adata = None
                redislogger.info(job_id, "AnnData object for MAGIC imputation is saved successfully")
                imputation_results['datasetId'] = datasetId
                create_pp_results(process_id, imputation_results)  # Insert pre-process results to database
                process_ids.append(process_id)
            else:
                adata = load_anndata(input)
                if adata is None:
                    detail = f"The file format is not supported: {input}"
                    upsert_jobs(
                        {
                            "job_id": job_id, 
                            "results": detail,
                            "completed_on": datetime.now(),
                            "status": "Failure"
                        }
                    )
                    raise CeleryTaskException(detail)
                elif 'MAGIC' not in adata.layers.keys(): 
                    try:
                        redislogger.info(job_id, "Start Magic imputation...")
                        counts = adata.X
                        data_magic = magic_impute(counts, genes)
                        adata.layers['MAGIC'] = data_magic
                        adata.write_h5ad(output, compression='gzip')

                        redislogger.info(job_id, "Computing PCA, neighborhood graph, tSNE, UMAP, and 3D UMAP")
                        adata, msg = run_dimension_reduction(adata, layer='MAGIC', n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=random_state)
                        if msg is not None: redislogger.warning(job_id, msg)

                        redislogger.info(job_id, "Clustering the neighborhood graph.")
                        adata = run_clustering(adata, layer='MAGIC', resolution=resolution, random_state=random_state)

                        redislogger.info(job_id, "Retrieving metadata and embeddings from AnnData object.")
                        imputation_results = get_metadata_from_anndata(adata, pp_stage, process_id, process, method, parameters, md5, adata_path=output, layer='MAGIC')
                        adata.write_h5ad(output, compression='gzip')
                        imputation_output.append({"MAGIC": output})
                        imputation_results["outputs"] = imputation_output
                        adata = None
                        redislogger.info(job_id, "AnnData object for MAGIC imputation is saved successfully")
                        process_ids.append(process_id)

                    except Exception as e:
                        detail = f"MAGIC imputation is failed: {e}"
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
                else: 
                    redislogger.warning(job_id, "'MAGIC' layer already exists.")
                    imputation_results = get_metadata_from_anndata(adata, pp_stage, process_id, process, method, parameters, md5, adata_path=output, layer='MAGIC')

            imputation_results['datasetId'] = datasetId
            create_pp_results(process_id, imputation_results)  # Insert pre-process results to database

        pp_results.append(imputation_results)
        

    # if "scGNN" in methods:
    #     if 'scGNN_imputed' not in adata.layers.keys(): 
    #         try:
    #             output = get_output_path(output, process_id, dataset, method='scGNN_imputation')
    #             redislogger.info(job_id, "AnnData object for scGNN imputation is saved successfully")          
    #         except Exception as e:
    #             redislogger.error(job_id, "scGNN imputation is failed.")
    #             if show_error: redislogger.error(job_id, f"scGNN imputation is failed: {e}")
    #     else: 
    #         redislogger.warning(job_id, "'scGNN_imputed' layer already exists.") 
    

    if "SAVER" in methods:
        method='SAVER'
        process_id = generate_process_id(md5, process, method, parameters)
        imputation_results = pp_result_exists(process_id)
        adata = None
        output = get_output_path(output, process_id, dataset, method='SAVER_imputation')
        csv_path = output.replace(".h5ad", ".csv")

        if imputation_results is not None:
            redislogger.info(job_id, "Found existing pre-process results in database, skip SAVER imputation.")
            process_ids.append(process_id)
            imputation_output = imputation_results["outputs"]
        else:
            if os.path.exists(output): # If output exist from the last run, then just pick up it.
                redislogger.info(job_id, "Output already exists, start from the last run.")
                adata = load_anndata(output)
                redislogger.info(job_id, "Computing PCA, neighborhood graph, tSNE, UMAP, and 3D UMAP")
                adata, msg = run_dimension_reduction(adata, layer='SAVER', n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=random_state)
                if msg is not None: redislogger.warning(job_id, msg)

                redislogger.info(job_id, "Clustering the neighborhood graph.")
                adata = run_clustering(adata, layer='SAVER', resolution=resolution, random_state=random_state)

                redislogger.info(job_id, "Retrieving metadata and embeddings from AnnData object.")
                imputation_results = get_metadata_from_anndata(adata, pp_stage, process_id, process, method, parameters, md5, adata_path=output, layer='SAVER')

                adata.write_h5ad(output, compression='gzip')
                imputation_output.append({"SAVER": output})
                imputation_output.append({"Report": report_path})
                imputation_results["outputs"] = imputation_output
                adata = None
                redislogger.info(job_id, "AnnData object for SAVER imputation is saved successfully")
                imputation_results['datasetId'] = datasetId
                create_pp_results(process_id, imputation_results)  # Insert pre-process results to database
                process_ids.append(process_id)
            else:
                adata, counts, csv_path = load_anndata_to_csv(input, csv_path)
                if adata is None:
                    detail = f"Layer {layer} does not exist in AnnData file: {input}"
                    raise CeleryTaskException(detail)
                elif 'SAVER' not in adata.layers.keys(): 
                    try:
                        # report_path = get_report_path(dataset, output, "SAVER")
                        report_path = output.replace(".h5ad", "_report.html")
                        
                        # Get the absolute path of the current file
                        current_file = os.path.abspath(__file__)

                        # Construct the relative path to the desired file
                        relative_path = os.path.join(os.path.dirname(current_file), 'imputation', 'SAVER.Rmd')

                        # Get the absolute path of the desired file
                        saver_path = os.path.abspath(relative_path)

                        redislogger.info(job_id, " Start SAVER imputation ...")

                        # saver_path = os.path.abspath("imputation/SAVER.Rmd")
                        s = subprocess.call([f"R -e \"rmarkdown::render('{saver_path}', params=list(dataset='{dataset}', input='{csv_path}', output='{output}', output_format='AnnData', ncores={ncores}), output_file='{report_path}')\""], shell = True)
                        # redislogger.info(job_id, str(s))

                        if os.path.exists(output):
                            adata = load_anndata(output)
                            redislogger.info(job_id, "Computing PCA, neighborhood graph, tSNE, UMAP, and 3D UMAP")
                            adata, msg = run_dimension_reduction(adata, layer='SAVER', n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=random_state)
                            if msg is not None: redislogger.warning(job_id, msg)

                            redislogger.info(job_id, "Clustering the neighborhood graph.")
                            adata = run_clustering(adata, layer='SAVER', resolution=resolution, random_state=random_state)

                            redislogger.info(job_id, "Retrieving metadata and embeddings from AnnData object.")
                            imputation_results = get_metadata_from_anndata(adata, pp_stage, process_id, process, method, parameters, md5, adata_path=output, layer='SAVER')

                            adata.write_h5ad(output, compression='gzip')
                            
                            imputation_output.append({"SAVER": output})
                            imputation_output.append({"Report": report_path})
                            imputation_results["outputs"] = imputation_output
                            adata = None
                            redislogger.info(job_id, "AnnData object for SAVER imputation is saved successfully")
                            process_ids.append(process_id)
                        else:
                            upsert_jobs(
                                {
                                    "job_id": job_id, 
                                    "results": "AnnData file does not exist due to the failure of SAVER Imputation.",
                                    "completed_on": datetime.now(),
                                    "status": "Failure"
                                }
                            )
                            raise ValueError("AnnData file does not exist due to the failure of SAVER Imputation.")
                
                    except Exception as e:
                        detail = f"SAVER imputation is failed: {e}"
                        upsert_jobs(
                            {
                                "job_id": job_id, 
                                "results": detail,
                                "completed_on": datetime.now(),
                                "status": "Failure"
                            }
                        )
                        raise CeleryTaskException(detail)
                else: 
                    redislogger.warning(job_id, "'SAVER' layer already exists.")
                    imputation_results = get_metadata_from_anndata(adata, pp_stage, process_id, process, method, parameters,  md5, adata_path=output, layer='SAVER')
                imputation_results['datasetId'] = datasetId
                create_pp_results(process_id, imputation_results)  # Insert pre-process results to database

        pp_results.append(imputation_results)
        
    process_ids = list(set(process_ids)) # De-duplicate process_ids

    results = {
        "output": imputation_output,
        "md5": md5,
        "process_ids": process_ids,
    }

    upsert_jobs(
        {
            "job_id": job_id, 
            "output": imputation_output,
            "results": results,
            "completed_on": datetime.now(),
            "status": "Success"
        }
    )

    return results

    
    

        
            

