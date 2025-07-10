import os
import subprocess
from tools.formating.formating import *
from config.celery_utils import get_input_path, get_output
from utils.redislogger import *
from utils.unzip import unzip_file_if_compressed
from fastapi import HTTPException, status
from tools.integration.scvi import scvi_integrate
from tools.reduction.reduction import run_dimension_reduction, run_clustering
from utils.mongodb import generate_process_id, pp_result_exists, create_pp_results, upsert_jobs
from exceptions.custom_exceptions import CeleryTaskException
from datetime import datetime

import warnings
warnings.simplefilter("ignore", FutureWarning)
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

def run_integration(job_id, ids:dict, fig_path=None):
    pp_stage = "Corrected"
    md5 = []
    process_ids = []
    process = "Integration"
    datasetIds = ids['datasetIds']
    datasets = ids['dataset']
    inputs = ids['input']
    userID = ids['userID']
    output = ids['output']
    do_umap = ids['do_umap']
    do_cluster = ids['do_cluster']
    
    # output_format = ids['output_format']
    parameters = ids['integration_params']
    batch_key = parameters['batch_key']
    pseudo_replicates = parameters['pseudo_replicates']
    methods = parameters['methods']
    default_assay = parameters['default_assay']
    # reference = parameters['reference']
    dims = parameters['dims']
    npcs = parameters['npcs']
    resolution = parameters['resolution']
    integration_output = []
    adata_outputs = []

    upsert_jobs(
        {
            "job_id": job_id, 
            "created_by": userID,
            "status": "Processing"
        }
    )

    if methods is None:
        redislogger.warning(job_id, "No integration method is selected.")
        detail = 'No integration method is selected.'
        raise CeleryTaskException(detail)
    # output = get_output_path(datasets, input, method='integration')
    methods = [x.upper() for x in methods if isinstance(x,str)]
    # adata, counts, csv_path = LoadAnndata_to_csv(input, output, layer, show_error)

    # methods = list_py_to_r(methods)
    abs_inputList = []

    if inputs is not None:
        for input in inputs:
            if input is not None:
                input = unzip_file_if_compressed(job_id, input)
                md5 = md5 + get_md5(input)
                abs_inputList.append(input)
    else:
        raise CeleryTaskException("No input file is found.")

    if datasets is not None:
        dataset = datasets[0]
    datasets = list_to_string(datasets)
    input = list_to_string_default(abs_inputList)
    
    redislogger.info(job_id, f"Using Integration Parameters: {parameters}")

    # #Get the absolute path for the given input
    # input = get_input_path(input, userID)
    # Get the absolute path for the given output
    # output = get_output(output, userID, job_id)
    for method in methods:
        process_id = generate_process_id(md5, process, method, parameters)
        output = os.path.join(os.path.dirname(inputs[0]), 'integration', f"{method}_integration.h5seurat")
        if not os.path.exists(os.path.dirname(output)):
            os.makedirs(os.path.dirname(output))
        adata_path = output.replace(".h5seurat", ".h5ad")
        report_path = output.replace(".h5seurat", ".html")

        # output = get_output_path(output, 'Integration', dataset=dataset, method=method, format='Seurat')
        # adata_path = get_output_path(output, 'Integration', dataset=dataset, method=method, format='AnnData')

        integration_results = pp_result_exists(process_id)

        if integration_results is not None:
            redislogger.info(job_id, "Found existing pre-process results in database, skip Integration.")
            integration_output = integration_results['outputs']
            adata_outputs.append(integration_results['adata_path'])
            process_ids.append(process_id)
        else:
            try:
                if ("HARMONY" in methods or "SCVI" in methods) and batch_key is not None :
                    adata = None
                    if len(inputs) > 1:
                        adatas = [load_anndata(input) for input in inputs]
                        for input in inputs:
                            ad = load_anndata(input)
                            if batch_key is None or batch_key.strip() == '':
                                ad.obs['batch'] = os.path.basename(input).split('.')[0] # If bacth_key is empty, use filename as batch_key                           
                        adata = sc.concat(adatas, join='outer')
                        if batch_key is None or batch_key.strip() == '':
                            batch_key = 'batch'
                    elif len(inputs) == 1:
                        adata = load_anndata(inputs[0])

                    # Check if X is normalized
                    if adata is not None:
                        if not is_normalized(adata.X) and check_nonnegative_integers(adata.X):
                            adata.layers['raw_counts'] = adata.X.copy() # Keep a copy of the raw counts
                    else:
                        raise CeleryTaskException(f"{method} integration is failed: AnnData is None.")
                    
                    # Pseudo replicates
                    if pseudo_replicates > 1:
                        adata = create_pseudo_replicates(adata, batch_key, pseudo_replicates)
      
                    if "HARMONY" in methods and adata is not None:
                        redislogger.info(job_id, "Start Harmony integration...")
                        import scanpy.external as sce
                        sc.pp.normalize_total(adata)
                        sc.pp.log1p(adata)
                        sc.pp.highly_variable_genes(adata, batch_key = batch_key, subset=False)
                        sc.pp.scale(adata)
                        sc.pp.pca(adata, use_highly_variable=True) #True since we didnt subset
                        sce.pp.harmony_integrate(adata, key = batch_key)
                        sc.pp.neighbors(adata, use_rep = "X_pca_harmony")
                        if do_umap:
                            redislogger.info(job_id, "Computing PCA, neighborhood graph, tSNE, UMAP, and 3D UMAP")
                            adata, msg = run_dimension_reduction(adata, n_neighbors=dims, n_pcs=npcs, use_rep="X_pca_harmony", random_state=0)
                            if msg is not None: redislogger.warning(job_id, msg)
                        if do_cluster:
                            redislogger.info(job_id, "Clustering the neighborhood graph.")
                            adata = run_clustering(adata, resolution=resolution, use_rep="X_pca_harmony", random_state=0)

                        redislogger.info(job_id, "Retrieving metadata and embeddings from AnnData object.")
                        integration_results = get_metadata_from_anndata(adata, pp_stage, process_id, process, method, parameters, md5, adata_path=adata_path, scanpy_cluster=batch_key)
                        adata.write_h5ad(output, compression='gzip')

                        integration_output.append({f"{method}_AnnDate": adata_path})
                        integration_results['outputs'] = integration_output
                        adata_outputs.append(adata_path)
                        adata = None
                        redislogger.info(job_id, integration_results['info'])
                        integration_results['datasetIds'] = datasetIds
                        create_pp_results(process_id, integration_results)  # Insert pre-process results to database 
                        process_ids.append(process_id) 

                    if "SCVI" in methods and adata is not None:
                        redislogger.info(job_id, "Start scVI integration...")
                        scvi_path = get_scvi_path(adata_path, "batch_integration")
                        adata.X = adata.layers['raw_counts'].copy() # Restore the raw counts

                        adata = scvi_integrate(adata, batch_key=batch_key, model_path=scvi_path)

                        sc.pp.neighbors(adata, use_rep = "X_scVI")
                        if do_umap:
                            redislogger.info(job_id, "Computing PCA, neighborhood graph, tSNE, UMAP, and 3D UMAP")
                            adata, msg = run_dimension_reduction(adata, n_neighbors=dims, n_pcs=npcs, use_rep="X_scVI", random_state=0)
                            if msg is not None: redislogger.warning(job_id, msg)
                        if do_cluster:
                            redislogger.info(job_id, "Clustering the neighborhood graph.")
                            adata = run_clustering(adata, resolution=resolution, use_rep="X_scVI", random_state=0)

                        redislogger.info(job_id, "Retrieving metadata and embeddings from AnnData object.")
                        integration_results = get_metadata_from_anndata(adata, pp_stage, process_id, process, method, parameters, md5, adata_path=adata_path, scanpy_cluster=batch_key)
                        adata.write_h5ad(output, compression='gzip')

                        integration_output.append({f"{method}_AnnDate": adata_path})
                        integration_results['outputs'] = integration_output
                        adata_outputs.append(adata_path)
                        adata = None
                        redislogger.info(job_id, integration_results['info'])
                        integration_results['datasetIds'] = datasetIds
                        create_pp_results(process_id, integration_results)  # Insert pre-process results to database 
                        process_ids.append(process_id)

                else:
                    redislogger.info(job_id, f"Start {method} integration...")
                    # report_path = get_report_path(dataset, output, "integration")
                    # Get the absolute path of the current file
                    current_file = os.path.abspath(__file__)
                    # Construct the relative path to the desired file
                    relative_path = os.path.join(os.path.dirname(current_file), 'integration', 'integration.Rmd')
                    # Get the absolute path of the desired file
                    rmd_path = os.path.abspath(relative_path)
                    # s = subprocess.call([f"R -e \"rmarkdown::render('{rmd_path}', params=list(unique_id='{job_id}', datasets='{datasets}', inputs='{input}', output_folder='{output}', adata_path='{adata_path}', methods='{methods}', dims='{dims}', npcs='{npcs}', default_assay='{default_assay}', reference='{reference}'), output_file='{report_path}')\""], shell = True)
                    s = subprocess.call([f"R -e \"rmarkdown::render('{rmd_path}', params=list(unique_id='{job_id}', datasets='{datasets}', inputs='{input}', output_folder='{output}', adata_path='{adata_path}', methods='{method}', dims={dims}, npcs={npcs}, resolution={resolution}, default_assay='{default_assay}'), output_file='{report_path}')\""], shell = True)
                    # redislogger.info(job_id, str(s))
                    # print(f"R -e \"rmarkdown::render('{rmd_path}', params=list(unique_id='{job_id}', datasets='{datasets}', inputs='{input}', output_folder='{output}', adata_path='{adata_path}', methods='{method}', dims={dims}, npcs={npcs}, default_assay='{default_assay}'), output_file='{report_path}')\"")

                    if os.path.exists(adata_path):
                        redislogger.info(job_id, "Adding 2D & 3D UMAP to AnnData object.")
                        adata = load_anndata(adata_path)
                        sc.pp.neighbors(adata, n_neighbors=dims, n_pcs=npcs, random_state=0)
                        adata = sc.tl.umap(adata, random_state=0, 
                                        init_pos="spectral", n_components=2, 
                                        copy=True, maxiter=None)
                        adata_3D = sc.tl.umap(adata, random_state=0, 
                                        init_pos="spectral", n_components=3, 
                                        copy=True, maxiter=None)
                        adata.obsm["X_umap_3D"] = adata_3D.obsm["X_umap"]

                        # Pseudo replicates
                        if pseudo_replicates > 1:
                            adata = create_pseudo_replicates(adata, batch_key, pseudo_replicates)

                        # Converrt dense martrix to sparse matrix
                        if isinstance(adata.X, np.ndarray):
                            adata.X = csr_matrix(adata.X)

                        adata.write_h5ad(adata_path, compression='gzip')
                        adata_3D = None
                    else:
                        upsert_jobs(
                            {
                                "job_id": job_id, 
                                "results": "AnnData file does not exist due to the failure of Integration.",
                                "completed_on": datetime.now(),
                                "status": "Failure"
                            }
                        )
                        raise ValueError("AnnData file does not exist due to the failure of Integration.")
                
                    redislogger.info(job_id, "Retrieving metadata and embeddings from AnnData object.")
                    integration_results = get_metadata_from_anndata(adata, pp_stage, process_id, process, method, parameters, md5, adata_path=adata_path, seurat_path=output, scanpy_cluster=batch_key)
                    # integration_output.append({method: {'adata_path': adata_path, 'seurat_path': output}})
                    integration_output.append({f"{method}_AnnDate": adata_path})
                    integration_output.append({f"{method}_Seurat": output})
                    integration_output.append({f"{method}_Report": report_path})
                    integration_results['outputs'] = integration_output
                    adata_outputs.append(adata_path)
                    adata = None
                    redislogger.info(job_id, integration_results['info'])
                    integration_results['datasetIds'] = datasetIds
                    create_pp_results(process_id, integration_results)  # Insert pre-process results to database 
                    process_ids.append(process_id) 

            except Exception as e:
                upsert_jobs(
                    {
                        "job_id": job_id, 
                        "results": f"Integration is failed: {e}",
                        "completed_on": datetime.now(),
                        "status": "Failure"
                    }
                )
                redislogger.error(job_id, f"{method} integration is failed: {e}")
                raise CeleryTaskException(f"{method} integration is failed: {e}")

    results = {
        "output": integration_output,
        "default_assay": default_assay,
        "md5": md5,
        "adata_path": adata_outputs,
        "process_ids": process_ids
    }

    upsert_jobs(
        {
            "job_id": job_id, 
            "output": integration_output,
            "datasetIds": datasetIds,
            "process_ids": process_ids,
            "results": results,
            "completed_on": datetime.now(),
            "status": "Success"
        }
    )

    return results
