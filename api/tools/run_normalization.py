import os
import subprocess
# import sys
# sys.path.append('..')
from tools.formating.formating import *
from config.celery_utils import get_input_path, get_output
from utils.redislogger import *
from tools.reduction.reduction import run_dimension_reduction, run_clustering
from utils.mongodb import generate_process_id, pp_result_exists, create_pp_results, upsert_jobs
from utils.unzip import unzip_file_if_compressed
from fastapi import HTTPException, status
from exceptions.custom_exceptions import CeleryTaskException
from datetime import datetime


def run_normalization(job_id, ds:dict, random_state=0, show_error=True):

    pp_results = []
    process_ids = []
    normalization_output = []
    pp_stage = "Normalized"
    process = "Normalization"
    dataset = ds['dataset']
    input = ds['input']
    userID = ds['userID']
    output = ds['output']
    species = ds['species']
    idtype = ds['idtype']
    cluster_label = ds['cluster_label']
    output_format = ds['output_format']
    parameters = ds['normalization_params']
    methods = parameters['methods']
    default_assay = parameters['assay']
    n_neighbors = parameters['n_neighbors']
    n_pcs = parameters['n_pcs']
    resolution = parameters['resolution']
    datasetId = ds['dataset_id']
    status = 'Successful'
    failed_methods = []

    upsert_jobs(
        {
            "job_id": job_id, 
            "created_by": userID,
            "status": "Processing"
        }
    )

    if len(methods) <1:
        redislogger.error(job_id, "No normalization method is selected.")
        detail = 'No normalization method is selected.'
        raise CeleryTaskException(detail)
    print(f"Slected methods: {methods}")
    # Get the absolute path for the given input
    # input = get_input_path(input, userID)
    input = unzip_file_if_compressed(job_id, ds['input'])
    md5 = get_md5(input)
    # Get the absolute path for the given output
    # output = get_output(output, userID, job_id)

    # methods = [x.upper() for x in methods if isinstance(x,str)]
    seurat_path = get_output_path(output, "normalized", dataset, method='normalization', format='Seurat')
    adata_path = get_output_path(output, "normalized", dataset, method='normalization', format='AnnData')
    adata_sct_path = adata_path.replace(".h5ad", "_SCT.h5ad")
    
    # methods = list_py_to_r(methods)
    if os.path.exists(seurat_path): # If seurat_path exist from the last run, then just pick up it.
        input = seurat_path
        redislogger.info(job_id, "Output already exists, start from the last run.")

    # Check if there is existing pre-process results
    methods_to_remove = []
    for method in methods:
        process_id = generate_process_id(md5, process, method, parameters)
        normalization_results = pp_result_exists(process_id)
        print(f"{method}: {process_id}")

        if normalization_results is not None:
            redislogger.info(job_id, f"Found existing pre-process results in database, skip {method} normalization.")
            pp_results.append(normalization_results)
            process_ids.append(process_id)
            methods_to_remove.append(method)

    print(f"Methods to remove: {methods_to_remove}.")
    if len(methods_to_remove) > 1:
        for method in methods_to_remove:
            methods.remove(method) # Remove method from methods list
    print(f"Remaining methods: {methods}.")

    if len(methods) > 0:
        try:
            report_path = get_report_path(dataset, output, "normalization")

            # Get the absolute path of the current file
            current_file = os.path.abspath(__file__)

            # Construct the relative path to the desired file
            relative_path = os.path.join(os.path.dirname(current_file), 'normalization', 'normalization.Rmd')

            # Get the absolute path of the desired file
            rmd_path = os.path.abspath(relative_path)

            # rmd_path = os.path.abspath("normalization/normalization.Rmd")
            s = subprocess.call([f"R -e \"rmarkdown::render('{rmd_path}', params=list(unique_id='{job_id}', dataset='{dataset}', input='{input}', output='{seurat_path}', adata_path='{adata_path}', output_format='{output_format}', methods='{list_to_string(methods)}', default_assay='{default_assay}', species='{species}', idtype='{idtype}'), output_file='{report_path}')\""], shell = True)
            # redislogger.info(job_id, str(s))

            if os.path.exists(adata_path):
                adata = load_anndata(adata_path)
                for layer in adata.layers.keys():
                    if layer in methods_to_remove: # Skip existing layers
                        continue

                    method = layer
                    process_id = generate_process_id(md5, process, method, parameters)
                    if method != "SCTransform":
                        try:
                            redislogger.info(job_id, f"Computing PCA, neighborhood graph, tSNE, UMAP, and 3D UMAP for layer {layer}.")
                            adata, msg = run_dimension_reduction(adata, layer=layer, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=random_state)
                            if msg is not None: redislogger.warning(job_id, msg)

                            redislogger.info(job_id, f"Clustering the neighborhood graph for layer {layer}.")
                            adata = run_clustering(adata, layer=layer, resolution=resolution, random_state=random_state)

                            redislogger.info(job_id, f"Retrieving metadata and embeddings from AnnData layer {layer}.")
                            normalization_results = get_metadata_from_anndata(adata, pp_stage, process_id, process, method, parameters, md5, layer=layer, adata_path=adata_path, seurat_path=output, cluster_label=cluster_label)
                            pp_results.append(normalization_results)
                            process_ids.append(process_id)
                            normalization_results['datasetId'] = datasetId
                            create_pp_results(process_id, normalization_results)  # Insert pre-process results to database
                        except Exception as e:
                            redislogger.error(job_id, f"UMAP or clustering is failed for {layer}: {e}")
                            failed_methods.append(f"UMAP or clustering is failed for {method}: {e}")
                    else:
                        if os.path.exists(adata_sct_path):
                            adata_sct = load_anndata(adata_sct_path)
                            method = "SCTransform"
                            process_id = generate_process_id(md5, process, method, parameters)
                            try:
                                redislogger.info(job_id, f"Computing PCA, neighborhood graph, tSNE, UMAP, and 3D UMAP for {method} normalization..")
                                adata_sct, msg = run_dimension_reduction(adata_sct, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=random_state)
                                if msg is not None: redislogger.warning(job_id, msg)

                                redislogger.info(job_id, f"Clustering the neighborhood graph for {method} normalization.")
                                adata_sct = run_clustering(adata_sct, resolution=resolution, random_state=random_state)

                                redislogger.info(job_id, f"Retrieving metadata and embeddings from AnnData normalized by {method}.")
                                normalization_results = get_metadata_from_anndata(adata_sct, pp_stage, process_id, process, method, parameters, md5, adata_path=adata_sct_path, seurat_path=output, cluster_label=cluster_label)
                                pp_results.append(normalization_results)
                                process_ids.append(process_id)
                                
                                adata.write_h5ad(adata_sct_path, compression='gzip')
                                normalization_results['datasetId'] = datasetId
                                create_pp_results(process_id, normalization_results)  # Insert pre-process results to database
                            except Exception as e:
                                redislogger.error(job_id, f"UMAP or clustering is failed for {layer}: {e}")
                                failed_methods.append(f"UMAP or clustering is failed for {method}: {e}")
                adata.write_h5ad(adata_path, compression='gzip')
  
            print(failed_methods)
            
        except Exception as e:
            # redislogger.error(job_id, "Normalization is failed.")
            detail = f"Normalization is failed: {e}"
            upsert_jobs(
                {
                    "job_id": job_id, 
                    "results": detail,
                    "completed_on": datetime.now(),
                    "status": "Failure"
                }
            )
            raise CeleryTaskException(detail)
        normalization_output.append({'AnnData_path': adata_path})
        normalization_output.append({'Seurat_path': output})
        normalization_output.append({'Report': report_path})
        if os.path.exists(adata_sct_path): normalization_output.append({'adata_sct_path': adata_sct_path})

    results = {
            "output": normalization_output,
            "default_assay": default_assay,
            "md5": md5,
            "process_ids": process_ids,
            "failed_methods": failed_methods
        }
    
    upsert_jobs(
        {
            "job_id": job_id, 
            "datasetId": datasetId,
            "process_ids": process_ids,
            "output": normalization_output,
            "results": results,
            "completed_on": datetime.now(),
            "status": "Success"
        }
    )

    return results