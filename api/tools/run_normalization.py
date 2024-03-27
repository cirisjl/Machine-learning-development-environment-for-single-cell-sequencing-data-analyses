import os
import subprocess
# import sys
# sys.path.append('..')
from tools.formating.formating import *
from config.celery_utils import get_input_path, get_output
from utils.redislogger import *
from tools.reduction.reduction import run_dimension_reduction, run_clustering
from utils.mongodb import generate_process_id, pp_results_exists, create_pp_results
from utils.unzip import unzip_file_if_compressed


def run_normalization(task_id, ds:dict, random_state=0, show_error=True):
    results = []
    pp_results = []
    process_ids = []
    pp_stage = "Normalized"
    process = "Normalization"
    dataset = ds.dataset
    default_assay = ds.assay
    input = ds.input
    userID = ds.userID
    output = ds.output
    methods = ds.methods
    species = ds.species
    idtype = ds.idtype
    parameters = ds.normalization_params
    n_neighbors = parameters.n_neighbors
    n_pcs = parameters.n_pcs
    resolution = parameters.resolution
    if methods is None:
        redislogger.warning(task_id, "No normalization method is selected.")
        return None
    
    # Get the absolute path for the given input
    # input = get_input_path(input, userID)
    input = unzip_file_if_compressed(input)
    md5 = get_md5(input)
    # Get the absolute path for the given output
    output = get_output(output, userID, task_id)

    # methods = [x.upper() for x in methods if isinstance(x,str)]
    output = get_output_path(dataset, output, method='normalization', format='Seurat')
    adata_path = get_output_path(dataset, output, method='normalization', format='AnnData')
    # methods = list_py_to_r(methods)

    # Check if there is existing pre-process results
    for method in methods:
        process_id = generate_process_id(md5, process, method, parameters)
        normalization_results = pp_results_exists(process_id)

        if normalization_results is not None:
            redislogger.info(task_id, f"Found existing pre-process results in database, skip {method} normalization.")
            pp_results.append(normalization_results)
            process_ids.append(process_id)
            methods.remove(method) # Remove method from methods list

    if len(methods) > 0:
        methods = list_to_string(methods)
        try:
            report_path = get_report_path(dataset, output, "normalization")

            # Get the absolute path of the current file
            current_file = os.path.abspath(__file__)

            # Construct the relative path to the desired file
            relative_path = os.path.join(os.path.dirname(current_file), 'normalization', 'normalization.Rmd')

            # Get the absolute path of the desired file
            rmd_path = os.path.abspath(relative_path)

            # rmd_path = os.path.abspath("normalization/normalization.Rmd")
            s = subprocess.call(["R -e \"rmarkdown::render('" + rmd_path + "', params=list(unique_id='" + task_id + "', dataset='" + str(dataset) + "', input='" + input + "', output='" + output + "', adata_path='" + adata_path + "', output_format='" + output_format + "', methods='" + methods + "', default_assay='" + default_assay + "', species='" + str(species) + "', idtype='" + str(idtype) + "'), output_file='" + report_path + "')\""], shell = True)
            redislogger.info(task_id, s)

            if os.path.exists(adata_path):
                adata = load_anndata(adata_path)
                for layer in adata.layers.keys():
                    method=layer
                    process_id = generate_process_id(md5, process, method, parameters)

                    redislogger.info(task_id, "Computing PCA, neighborhood graph, tSNE, UMAP, and 3D UMAP")
                    adata, msg = run_dimension_reduction(adata, layer=layer, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=random_state)
                    if msg is not None: redislogger.warning(task_id, msg)

                    redislogger.info(task_id, "Clustering the neighborhood graph.")
                    adata = run_clustering(adata, layer=layer, resolution=resolution, random_state=random_state)

                    redislogger.info(task_id, "Retrieving metadata and embeddings from AnnData object.")
                    normalization_results = get_metadata_from_anndata(adata, pp_stage, process_id, process, method, parameters, adata_path, seurat_path=output)
                    pp_results.append(normalization_results)
                    process_ids.append(process_id)
                    create_pp_results(normalization_results)  # Insert pre-process results to database  
        except Exception as e:
            redislogger.error(task_id, "Normalization is failed.")
            if show_error: redislogger.error(task_id, f"Normalization is failed: {e}")

    results.append({
            "task_id": task_id, 
            "inputfile": input,
            "default_assay": default_assay,
            "md5": md5,
            "process_id": process_ids,
            "pp_results": pp_results,
            "message": "Normalization completed successfully."
        })  

    return results