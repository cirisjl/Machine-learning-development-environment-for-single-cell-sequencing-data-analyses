import os
import subprocess
from tools.formating.formating import *
from config.celery_utils import get_input_path, get_output
from utils.redislogger import *
from utils.unzip import unzip_file_if_compressed
from fastapi import HTTPException, status
from tools.reduction.reduction import run_dimension_reduction, run_clustering
from utils.mongodb import generate_process_id, pp_results_exists, create_pp_results, upsert_task_results


def run_integration(task_id, ids:dict):
    pp_stage = "Corrected"
    md5 = []
    process_ids = []
    process = "Integration"
    datasetIds = ids['datasetIds']
    datasets = ids['datasets']
    inputs = ids['inputs']
    userID = ids['userID']
    output = ids['output']
    methods = ids['methods']
    output_format = ids['output_format']
    parameters = ids['params']
    species = parameters['species']
    default_assay = parameters['default_assay']
    genes = parameters['genes']
    reference = parameters['']
    dims = parameters['dims']
    npcs = parameters['npcs']
    integration_output = []

    if methods is None:
        redislogger.warning(task_id, "No integration method is selected.")
        raise HTTPException(
            status_code = status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail = 'No integration method is selected.')
    # output = get_output_path(datasets, input, method='integration')
    # methods = [x.upper() for x in methods if isinstance(x,str)]
    # adata, counts, csv_path = LoadAnndata_to_csv(input, output, layer, show_error)

    # methods = list_py_to_r(methods)
    abs_inputList = []

    if inputs is not None:
        for input in inputs:
            if input is not None:
                input = unzip_file_if_compressed(task_id, input)
                md5 = md5 + get_md5(input)
                abs_inputList.append(input)

    if datasets is not None:
        dataset = datasets[0]
    datasets = list_to_string(datasets)
    methods = list_to_string(methods)
    input = list_to_string_default(abs_inputList)
    

    # #Get the absolute path for the given input
    # input = get_input_path(input, userID)
    #Get the absolute path for the given output
    # output = get_output(output, userID, task_id)
    for method in methods:
        process_id = generate_process_id(md5, process, method, parameters)
        output = get_output_path(output, 'Integration', dataset=dataset, method=method, format='Seurat')
        adata_path = get_output_path(output, 'Integration', dataset=dataset, method=method, format='AnnData')

        integration_results = pp_results_exists(process_id)

        if integration_results is not None:
            redislogger.info(task_id, "Found existing pre-process results in database, skip Integration.")
            integration_output.append({method: integration_results[method]})
        else:
            try:
                report_path = get_report_path(dataset, output, "integration")
                # Get the absolute path of the current file
                current_file = os.path.abspath(__file__)
                # Construct the relative path to the desired file
                relative_path = os.path.join(os.path.dirname(current_file), 'integration', 'integration.Rmd')
                # Get the absolute path of the desired file
                rmd_path = os.path.abspath(relative_path)
                s = subprocess.call(["R -e \"rmarkdown::render('" + rmd_path + "', params=list(unique_id='" + task_id + "', datasets='" + str(datasets) + "', inputs='" + str(input) + "', output_folder='" + output + "', adata_path='" + adata_path + "', output_format='" + output_format + "', methods='" + str(methods) + "', dims='" + str(dims) + "', npcs='" + str(npcs) + "', default_assay='" + default_assay + "', reference='" + str(reference) + "', genes='" + str(genes) + "'), output_file='" + report_path + "')\""], shell = True)
                redislogger.info(task_id, s)

                if os.path.exists(adata_path):
                    redislogger.info(task_id, "Adding 3D UMAP to AnnData object.")
                    adata = load_anndata(adata_path)
                    sc.pp.neighbors(adata, n_neighbors=dims, n_pcs=npcs, random_state=0)
                    adata_3D = sc.tl.umap(adata, random_state=0, 
                                    init_pos="spectral", n_components=3, 
                                    copy=True, maxiter=None)
                    adata.obsm["X_umap_3D"] = adata_3D.obsm["X_umap"]
                    adata.write_h5ad(adata_path, compression='gzip')
                    adata_3D = None
                else:
                    raise ValueError("AnnData file does not exist due to the failure of Bioconductor QC.")
                
                redislogger.info(task_id, "Retrieving metadata and embeddings from AnnData object.")
                integration_results = get_metadata_from_anndata(adata, pp_stage, process_id, process, method, ids, md5, adata_path=adata_path, seurat_path=output)
                integration_output.append({method: {'adata_path': adata_path, 'seurat_path': output}})
                adata = None
                redislogger.info(task_id, integration_results['info'])
                create_pp_results(process_id, integration_results)  # Insert pre-process results to database 
                process_ids.append(process_id) 

            except Exception as e:
                redislogger.error(task_id, f"Integration is failed: {e}")

    results = {
        "taskId": task_id, 
        "owner": userID,
        "inputfile": inputs,
        "output": integration_output,
        "default_assay": default_assay,
        "md5": md5,
        "process_ids": process_ids,
        # "pp_results": pp_results,
        "status": "Success"
    }

    upsert_task_results(results)

    return results
