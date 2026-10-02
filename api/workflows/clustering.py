import os
from tools.formating.formating import *
from tools.run_qc import run_qc
from tools.run_normalization import run_normalization
from tools.run_imputation import run_imputation
from utils.unzip import unzip_file_if_compressed
from fastapi import HTTPException, status
from utils.redislogger import *
from utils.mongodb import generate_workflow_id, upsert_jobs, upsert_workflows
from datetime import datetime
from exceptions.custom_exceptions import CeleryTaskException


def run_clustering(job_id, ds:dict, random_state=0):
    wf_results = {}
    process_ids = []
    userID = ds['userID']
    datasetId = ds['datasetId']
    dataset = ds['dataset']
    qc_params = ds['qc_params']
    normalization_params = ds['normalization_params']
    imputation_params = ds['imputation_params']
    fig_path = None
    reduction_params = ds['reduction_params']
    n_neighbors = reduction_params['n_neighbors']
    n_pcs = reduction_params['n_pcs']
    resolution = reduction_params['resolution']
    n_hvg = reduction_params['n_hvg']
    qc_params['n_neighbors'] = n_neighbors
    qc_params['n_pcs'] = n_pcs
    qc_params['resolution'] = resolution
    normalization_params['n_neighbors'] = n_neighbors
    normalization_params['n_pcs'] = n_pcs
    normalization_params['resolution'] = resolution
    imputation_params['n_neighbors'] = n_neighbors
    imputation_params['n_pcs'] = n_pcs
    imputation_params['resolution'] = resolution
    ds['n_hvg'] = n_hvg

    # Initialize methodMap
    methodMap = {}

    # Extract methods for Quality Control
    if qc_params and "methods" in qc_params and qc_params["methods"]:
        methodMap["Quality Control"] = ", ".join(qc_params["methods"])

    # Extract methods for Normalization
    if normalization_params and "methods" in normalization_params and normalization_params["methods"]:
        methodMap["Normalization"] = ", ".join(normalization_params["methods"])

    # Extract methods for Imputation
    if imputation_params and "methods" in imputation_params and imputation_params["methods"]:
        methodMap["Imputation"] = ", ".join(imputation_params["methods"])
        
    output = None
    description = "Clustering Workflow" 
    input = unzip_file_if_compressed(job_id, ds['input'])
    md5 = get_md5(input)

    # Create folders for output figures
    workflow_id = generate_workflow_id(md5, "clustering", ds)
    fig_path = os.path.join(os.path.dirname(ds['input']), 'workflow', workflow_id)
    if not os.path.exists(fig_path):
        os.makedirs(fig_path, exist_ok=True)
    wf_results['figures'] = fig_path

    # for method i want methodMap key value pairs.

    if datasetId is not None:
        description = f"Clustering Workflow for {datasetId}"
    elif dataset is not None:
        description = f"Clustering Workflow for {dataset}"

    # wf_results = pp_result_exists(process_id)
    upsert_jobs(
        {
            "job_id": job_id, 
            "created_by": userID,
            "Description": description,
            # "Method": str(methodMap).replace("'", "").replace("{", "").replace("}", ""),
            "Method": methodMap,
            "datasetURL": input,
            "datasetId": datasetId,
            "Process": "Clustering",
            "Category": 'workflow',
            "Created on": datetime.now(),
            "Status": "Processing"
        }
    )

    try:
        process_id = generate_process_id(md5, process, method, parameters)
        # methods = [x.upper() for x in methods if isinstance(x,str)]
        seurat_path = get_output_path(output, process_id, dataset, method=method, format='Seurat')
        adata_path = get_output_path(output, process_id, dataset, method=method, format='AnnData')
        # adata_sct_path = adata_path.replace(".h5ad", "_SCT.h5ad")
        report_path = adata_path.replace(".h5ad", "_report.html")
        
        # methods = list_py_to_r(methods)
        # if os.path.exists(seurat_path): # If seurat_path exist from the last run, then just pick up it.
        #     input = seurat_path
        #     redislogger.info(job_id, "Output already exists, start from the last run.")
        # report_path = get_report_path(dataset, output, "normalization")
        
        # Get the absolute path of the current file
        current_file = os.path.abspath(__file__)

        # Construct the relative path to the desired file
        relative_path = os.path.join(os.path.dirname(current_file), 'normalization', 'normalization.Rmd')

        # Get the absolute path of the desired file
        rmd_path = os.path.abspath(relative_path)

        # rmd_path = os.path.abspath("normalization/normalization.Rmd")
        # s = subprocess.call([f"R -e \"rmarkdown::render('{rmd_path}', params=list(unique_id='{job_id}', dataset='{dataset}', input='{input}', output='{seurat_path}', adata_path='{adata_path}', output_format='{output_format}', methods='{list_to_string(methods)}', default_assay='{default_assay}', species='{species}', idtype='{idtype}'), output_file='{report_path}')\""], shell = True)
        s = subprocess.call([f"R -e \"rmarkdown::render('{rmd_path}', params=list(unique_id='{job_id}', dataset='{dataset}', input='{input}', output='{seurat_path}', adata_path='{adata_path}', output_format='{output_format}', methods='{method}', default_assay='{default_assay}', species='{species}', idtype='{idtype}'), output_file='{report_path}')\""], shell = True)
        # redislogger.info(job_id, str(s))

        if os.path.exists(adata_path):
            adata = load_anndata(adata_path)
            # print(adata_path)
            # print(adata)
            # if method != "SCT" and method != "SCT V2":
            try:
                if do_umap:
                    redislogger.info(job_id, f"Computing PCA, neighborhood graph, tSNE, UMAP, and 3D UMAP for layer {method}.")
                    adata, msg = run_dimension_reduction(adata, layer=method, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=random_state)
                    if msg is not None: redislogger.warning(job_id, msg)
                if do_cluster:
                    redislogger.info(job_id, f"Clustering the neighborhood graph for layer {method}.")
                    adata = run_clustering(adata, layer=method, resolution=resolution, random_state=random_state)
                if fig_path is not None:
                    plot_embedding(adata, layer=method, fig_path=fig_path, title=method + " Normalization")
                
                # # Converrt dense martrix to sparse matrix
                # if isinstance(adata.X, np.ndarray):
                #     adata.X = csr_matrix(adata.X)
                # adata.write_h5ad(adata_path, compression='gzip')
                output_path, zarr_output = save_anndata(adata, adata_path, zarr=True, n_hvg=n_hvg, layer=method, obs_cols=obs_cols)

                redislogger.info(job_id, f"Retrieving metadata and embeddings from AnnData layer {method}.")
                normalization_results = get_metadata_from_anndata(adata, pp_stage, process_id, process, method, parameters, md5, layer=method, adata_path=adata_path, seurat_path=output, cluster_label=cluster_label, cluster_colname=cluster_colname, zarr_path=zarr_output)
                
                tools_path = adata_path.replace('.h5ad', '.txt')
                if os.path.exists(tools_path):
                    tools = get_r_tools(tools_path)
                    normalization_results['tools'] = tools
                if os.path.exists(adata_path): normalization_output.append({'AnnData': adata_path})
                if os.path.exists(seurat_path): normalization_output.append({'Seurat': seurat_path})
                if os.path.exists(report_path): normalization_output.append({'Report': report_path})

                # Add preset questions for tissue and species
                if organ_part is not None and organ_part != "" and species is not None and species != "":
                    if cluster_colname is None or cluster_colname.strip() == "":
                        cluster_colname = f"{method}_leiden"
                    preset_questions = create_annotation_prompt(adata, tissue=organ_part, species=species, layer=method, method="t-test", groupby=cluster_colname, top=n_hvg, task="Normalization")
                    normalization_results['preset_questions'] = preset_questions

                normalization_results['outputs'] = normalization_output
                adata = None
                
                # pp_results.append(normalization_results)
                process_ids.append(process_id)
                normalization_results['datasetId'] = datasetId
                normalization_results['created_by'] = userID
                create_pp_results(process_id, normalization_results)  # Insert pre-process results to database
                adata_paths.append(adata_path)
            except Exception as e:
                redislogger.error(job_id, f"UMAP or clustering is failed for {method}: {e}")
                failed_methods.append(f"UMAP or clustering is failed for {method}: {e}")
        else:
            redislogger.error(job_id, f'{method} normalization is failed, no output is created.')
            if not wf:
                raise Exception(f'{method} is failed, no output is created.')

        # if os.path.exists(adata_sct_path):
        #     adata_sct = load_anndata(adata_sct_path)
        #     method = "SCTransform"
        #     process_id = generate_process_id(md5, process, method, parameters)
        #     try:
        #         redislogger.info(job_id, f"Computing PCA, neighborhood graph, tSNE, UMAP, and 3D UMAP for {method} normalization..")
        #         adata_sct, msg = run_dimension_reduction(adata_sct, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=random_state)
        #         if msg is not None: redislogger.warning(job_id, msg)

        #         redislogger.info(job_id, f"Clustering the neighborhood graph for {method} normalization.")
        #         adata_sct = run_clustering(adata_sct, resolution=resolution, random_state=random_state, fig_path=fig_path)

        #         redislogger.info(job_id, f"Retrieving metadata and embeddings from AnnData normalized by {method}.")
        #         normalization_results = get_metadata_from_anndata(adata_sct, pp_stage, process_id, process, method, parameters, md5, adata_path=adata_sct_path, seurat_path=output, cluster_label=cluster_label)
        #         if os.path.exists(adata_sct_path): normalization_output.append({'AnnData_sct': adata_sct_path})
        #         normalization_results['outputs'] = normalization_output
                
        #         # pp_results.append(normalization_results)
        #         process_ids.append(process_id)
                
        #         adata_sct.write_h5ad(adata_sct_path, compression='gzip')
        #         normalization_results['datasetId'] = datasetId
        #         create_pp_results(process_id, normalization_results)  # Insert pre-process results to database
        #     except Exception as e:
        #         redislogger.error(job_id, f"UMAP or clustering is failed for SCTransform: {e}")
        #         failed_methods.append(f"UMAP or clustering is failed for SCTransform: {e}")

        print(failed_methods)
        
    except Exception as e:
        # redislogger.error(job_id, "Normalization is failed.")
        detail = f"Normalization is failed: {e}"
        upsert_jobs(
            {
                "job_id": job_id, 
                "results": detail,
                "Completed on": datetime.now(),
                "Status": "Failure"
            }
        )
        redislogger.error(job_id, detail)
        if not wf:
            raise CeleryTaskException(detail)

    try:
        qc_results = run_qc(job_id, ds, fig_path=fig_path, wf=True)
        if qc_results is not None:
            ds['input'] = qc_results['adata_path']
            print(qc_results['adata_path'])        
            process_ids.extend(qc_results["process_ids"])
            wf_results['QC'] = qc_results["process_ids"]
            wf_results['QC_output'] = qc_results['output']
            if normalization_params["methods"] is not None:
                normalization_results = run_normalization(job_id, ds, fig_path=fig_path, wf=True)
                wf_results['normalization'] = normalization_results["process_ids"]
                process_ids.extend(normalization_results["process_ids"])
                wf_results['normalization_output'] = normalization_results['output']
                output = normalization_results['output']
            elif imputation_params["methods"] is not None:
                imputation_results = run_imputation(job_id, ds, fig_path=fig_path, wf=True)
                wf_results['imputation'] = imputation_results["process_ids"]
                process_ids.extend(imputation_results["process_ids"])
                wf_results['imputation_output'] = imputation_results['output']
                output = imputation_results['output']
            # upsert_workflows(workflow_id, wf_results)

        results = {
            "output": output,
            # "workflow_id": workflow_id,
            "md5": md5,
            "wf_results": wf_results,
            # "figures":fig_path, 
            "process_ids": process_ids
        }
        
        upsert_jobs(
            {
                "job_id": job_id, 
                "output": output,
                "process_ids": process_ids,
                # "workflow_id": workflow_id,
                "results": results,
                # "figures": fig_path, 
                "Completed on": datetime.now(),
                "Status": "Success"
            }
        )
    except Exception as e:
        detail = f"Clustering workflow is failed: {e}"
        upsert_jobs(
            {
                "job_id": job_id, 
                "results": {"error": detail},
                "Completed on": datetime.now(),
                "Status": "Failure"
            }
        )
        raise CeleryTaskException(detail)

    return results 