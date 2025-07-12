import os
import subprocess
import sys
from tools.formating.formating import *
from tools.annotation.celltypist import run_celltypist
from tools.annotation.scvi import scvi_transfer
from config.celery_utils import get_input_path, get_output
from utils.redislogger import *
from tools.reduction.reduction import run_dimension_reduction, run_clustering
from utils.mongodb import generate_process_id, pp_result_exists, create_pp_results, upsert_jobs
from utils.unzip import unzip_file_if_compressed
from fastapi import HTTPException, status
from exceptions.custom_exceptions import CeleryTaskException
from datetime import datetime
    

def run_annotation(job_id, ds:dict, fig_path=None, show_error=True, random_state=0):
    pp_results = []
    process_ids = []
    annotation_output = []
    pp_stage = "Annotation"
    process = "Annotation"
    dataset = ds['dataset']
    species = ds['species'].lower()
    input = ds['input']
    user_refs = ds['user_refs']
    userID = ds['userID']
    output = ds['output']
    datasetId = ds['datasetId']
    do_umap = ds['do_umap']
    do_cluster = ds['do_cluster']
    parameters = ds['annotation_params']
    layer = None
    if parameters['layer'] is not None and parameters['layer'].strip != "":
        layer = parameters['layer']
    assay = parameters['assay']
    methods = parameters['methods']
    celltypist_model = parameters['celltypist_model']
    SingleR_ref = parameters['SingleR_ref']
    user_label = parameters['user_label']
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
        redislogger.error(job_id, "No annotation method is selected.")
        detail = 'No annotation method is selected.'
        raise CeleryTaskException(detail)
    
    input = unzip_file_if_compressed(job_id, ds['input'])
    md5 = get_md5(input)
    # process_id = generate_process_id(md5, process, methods, parameters)
    # output = get_output_path(output, process_id, dataset, method='annotation')
    adata_path = get_output_path(output, 'annotation', dataset)

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

    redislogger.info(job_id, f"Using Annotation Parameters: {parameters}")

    methods = [x.upper() for x in methods if isinstance(x,str)]
    for method in methods:
        process_id = generate_process_id(md5, process, method, parameters)
        # adata_path = get_output_path(output, process_id, dataset, method=method)
        annotation_results = pp_result_exists(process_id)

        if annotation_results is not None:
            redislogger.info(job_id, f"Found existing pre-process results in database, skip {method} annotation.")
            process_ids.append(process_id)
            annotation_output = annotation_results["outputs"]
        else:
            if method == "CELLTYPIST":
                try:
                    redislogger.info(job_id, "Start CellTypist annotation...")
                    adata = run_celltypist(adata, model_name=celltypist_model, refs = user_refs, labels = user_label, species = species)
                    redislogger.info(job_id, "CellTypist annotation has been added to AnnData.obs.")
                    
                    if do_umap:
                        redislogger.info(job_id, "Computing PCA, neighborhood graph, tSNE, UMAP, and 3D UMAP")
                        adata, msg = run_dimension_reduction(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=random_state)
                        if msg is not None: redislogger.warning(job_id, msg)
                    if do_cluster:
                        redislogger.info(job_id, "Clustering the neighborhood graph.")
                        adata = run_clustering(adata, resolution=resolution, random_state=random_state, fig_path=fig_path)

                    redislogger.info(job_id, "Retrieving metadata and embeddings from AnnData object.")
                    annotation_results = get_metadata_from_anndata(adata, pp_stage, process_id, process, method, parameters, md5, adata_path=adata_path)

                    # Converrt dense martrix to sparse matrix
                    if isinstance(adata.X, np.ndarray):
                        adata.X = csr_matrix(adata.X)
                    adata.write_h5ad(adata_path, compression='gzip')
                    annotation_output.append({"CellTypist": adata_path})
                    annotation_results["outputs"] = annotation_output
                    redislogger.info(job_id, "AnnData object for CellTypist annotation is saved successfully")
                    process_ids.append(process_id)

                    annotation_results['datasetId'] = datasetId
                    create_pp_results(process_id, annotation_results)  # Insert pre-process results to database
                    pp_results.append(annotation_results)

                except Exception as e:
                    detail = f"CellTypist annotation is failed: {e}"
                    upsert_jobs(
                        {
                            "job_id": job_id, 
                            "results": detail,
                            "completed_on": datetime.now(),
                            "status": "Failure"
                        }
                    )
                    # os.remove(output)
                    raise CeleryTaskException(detail)

            if method == "SCVI":
                try:
                    redislogger.info(job_id, "Start scVI annotation...")
                    adata = scvi_transfer(adata, refs = user_refs, labels = user_label)
                    redislogger.info(job_id, "scVI cell type transfer has been added to AnnData.obs.")

                    if do_umap:
                        redislogger.info(job_id, "Computing PCA, neighborhood graph, tSNE, UMAP, and 3D UMAP")
                        adata, msg = run_dimension_reduction(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=random_state)
                        if msg is not None: redislogger.warning(job_id, msg)
                    if do_cluster:
                        redislogger.info(job_id, "Clustering the neighborhood graph.")
                        adata = run_clustering(adata, resolution=resolution, random_state=random_state, fig_path=fig_path)

                    redislogger.info(job_id, "Retrieving metadata and embeddings from AnnData object.")
                    annotation_results = get_metadata_from_anndata(adata, pp_stage, process_id, process, method, parameters, md5, adata_path=adata_path)

                    # Converrt dense martrix to sparse matrix
                    if isinstance(adata.X, np.ndarray):
                        adata.X = csr_matrix(adata.X)
                    adata.write_h5ad(adata_path, compression='gzip')
                    annotation_output.append({"scVI": adata_path})
                    annotation_results["outputs"] = annotation_output
                    # adata = None
                    redislogger.info(job_id, "AnnData object for scVI annotation is saved successfully")
                    process_ids.append(process_id)
                    annotation_results['datasetId'] = datasetId
                    create_pp_results(process_id, annotation_results)  # Insert pre-process results to database
                    pp_results.append(annotation_results)

                except Exception as e:
                    detail = f"scVI annotation is failed: {e}"
                    upsert_jobs(
                        {
                            "job_id": job_id, 
                            "results": detail,
                            "completed_on": datetime.now(),
                            "status": "Failure"
                        }
                    )
                    # os.remove(adata_path)
                    raise CeleryTaskException(detail)

            if method == "SINGLER":
                if SingleR_ref is None and (len(refs) == 0 or labels is None):
                    raise CeleryTaskException(f"SingleR annotation is failed due to empty reference ({SingleR_ref}) and empty user reference ({user_refs}) or cell labels ({user_label}).")

                try:
                    # report_path = get_report_path(dataset, output, "SAVER")
                    report_path = adata_path.replace(".h5ad", "_report.html")
                    output_folder = os.path.dirname(adata_path)
                    
                    # Get the absolute path of the current file
                    current_file = os.path.abspath(__file__)

                    # Construct the relative path to the desired file
                    relative_path = os.path.join(os.path.dirname(current_file), 'annotation', 'singleR.Rmd')

                    # Get the absolute path of the desired file
                    singler_path = os.path.abspath(relative_path)

                    redislogger.info(job_id, " Start SingleR annotation ...")
                    if user_label is not None and len(user_refs) > 0:
                        s = subprocess.call([f"R -e \"rmarkdown::render('{singler_path}', params=list(unique_id='{job_id}', dataset='{dataset}', input='{input}', output_folder='{output_folder}', dims={n_neighbors}, npcs={n_pcs}, resolution={resolution}, species='{species}', default_assay='{assay}', reference='{SingleR_ref}', user_ref='{user_refs[0]}', user_label='{user_label}'), output_file='{report_path}')\""], shell = True)
                    else:
                        s = subprocess.call([f"R -e \"rmarkdown::render('{singler_path}', params=list(unique_id='{job_id}', dataset='{dataset}', input='{input}', output_folder='{output_folder}', dims={n_neighbors}, npcs={n_pcs}, resolution={resolution}, species='{species}', default_assay='{assay}', reference='{SingleR_ref}'), output_file='{report_path}')\""], shell = True)

                    csv_main = output_folder + "/results_main.csv"
                    csv_fine = output_folder + "/results_fine.csv"
                    csv_user = output_folder + "/results_user.csv"

                    if not (os.path.exists(csv_main) or os.path.exists(csv_fine) or os.path.exists(csv_user)):
                        upsert_jobs(
                            {
                                "job_id": job_id, 
                                "results": "SingleR annotation is failed.",
                                "completed_on": datetime.now(),
                                "status": "Failure"
                            }
                        )
                        # redislogger.warning(job_id, 'SingleR annotation is failed.')
                        raise CeleryTaskException('SingleR annotation is failed.')

                    if os.path.exists(csv_main):
                        df_main = pd.read_csv(csv_main, index_col=0)
                        adata.obs['SingleR_main'] = df_main['labels']
                        adata.obs['SingleR_main.pruned'] = df_main['pruned.labels']
                    
                    if os.path.exists(csv_fine):
                        df_fine = pd.read_csv(csv_fine, index_col=0)
                        adata.obs['SingleR_fine'] = df_fine['labels']
                        adata.obs['SingleR_fine.pruned'] = df_fine['pruned.labels']

                    if os.path.exists(csv_user):
                        df_user = pd.read_csv(csv_user, index_col=0)
                        adata.obs['SingleR_user_ref'] = df_user['labels']
                        adata.obs['SingleR_user_ref.pruned'] = df_user['pruned.labels']

                    if do_umap:
                        redislogger.info(job_id, "Computing PCA, neighborhood graph, tSNE, UMAP, and 3D UMAP")
                        adata, msg = run_dimension_reduction(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=random_state)
                        if msg is not None: redislogger.warning(job_id, msg)
                    if do_cluster:
                        redislogger.info(job_id, "Clustering the neighborhood graph.")
                        adata = run_clustering(adata, resolution=resolution, random_state=random_state, fig_path=fig_path)

                    redislogger.info(job_id, "Retrieving metadata and embeddings from AnnData object.")
                    annotation_results = get_metadata_from_anndata(adata, pp_stage, process_id, process, method, parameters, md5, adata_path=adata_path)
                    
                    # Converrt dense martrix to sparse matrix
                    if isinstance(adata.X, np.ndarray):
                        adata.X = csr_matrix(adata.X)
                    adata.write_h5ad(adata_path, compression='gzip')
                    
                    annotation_output.append({"SingleR": adata_path})
                    annotation_output.append({"Report": report_path})
                    annotation_results["outputs"] = annotation_output
                    redislogger.info(job_id, "AnnData object for SingleR annotation is saved successfully")
                    process_ids.append(process_id)
                    annotation_results['datasetId'] = datasetId
                    create_pp_results(process_id, annotation_results)  # Insert pre-process results to database
                    pp_results.append(annotation_results)
            
                except Exception as e:
                    detail = f"SingleR annotation is failed: {e}"
                    upsert_jobs(
                        {
                            "job_id": job_id, 
                            "results": detail,
                            "completed_on": datetime.now(),
                            "status": "Failure"
                        }
                    )
                    raise CeleryTaskException(detail)

        
    process_ids = list(set(process_ids)) # De-duplicate process_ids

    results = {
        "output": annotation_output,
        "adata_path": adata_path,
        "md5": md5,
        "process_ids": process_ids,
    }

    upsert_jobs(
        {
            "job_id": job_id, 
            "datasetId": datasetId,
            "process_ids": process_ids,
            "adata_path": adata_path,
            "output": annotation_output,
            "results": results,
            "completed_on": datetime.now(),
            "status": "Success"
        }
    )

    adata = None

    return results

    
    

        
            

