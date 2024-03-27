import os
import subprocess
import sys
from tools.formating.formating import *
from tools.imputation.MAGIC import magic_impute
from config.celery_utils import get_input_path, get_output
from utils.redislogger import *
from tools.utils.reduction import run_dimension_reduction, run_clustering
    

def run_imputation(task_id, dataset, input, userID, output, methods, layer=None, genes=None, ncores=12, show_error=True):
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
        redislogger.warning(task_id, "No imputation method is selected.")
        return None
    
    #Get the absolute path for the given input
    # input = get_input_path(input, userID)
    #Get the absolute path for the given output
    output = get_output(output, userID, task_id)
    methods = [x.upper() for x in methods if isinstance(x,str)]
    
    if "MAGIC" in methods:
        adata = load_anndata(input)
        if adata is None:
            redislogger.error(task_id, f"File format is not supported: {input}")
        elif 'MAGIC_imputed' not in adata.layers.keys(): 
            try:
                redislogger.info(task_id, "Start Magic imputation...")
                counts = adata.X
                data_magic = magic_impute(counts, genes)
                adata.layers['MAGIC_imputed'] = data_magic
                redislogger.info(task_id, "Computing PCA, neighborhood graph, tSNE, UMAP, and 3D UMAP")
                adata, msg = run_dimension_reduction(adata, n_neighbors=parameters.n_neighbors, n_pcs=parameters.n_pcs, random_state=random_state)
                if msg is not None: redislogger.warning(task_id, msg)

                redislogger.info(task_id, "Clustering the neighborhood graph.")
                adata = run_clustering(adata, resolution=parameters.resolution, random_state=random_state)

                redislogger.info(task_id, "Retrieving metadata and embeddings from AnnData object.")
                normalization_results = get_metadata_from_anndata(adata, pp_stage, process_id, process, method, parameters, adata_path, seurat_path=output)


                output = get_output_path(dataset, output, method='MAGIC_imputation')
                adata.write_h5ad(output, compression='gzip')
                redislogger.info(task_id, "AnnData object for MAGIC imputation is saved successfully")
            except Exception as e:
                redislogger.error(task_id, "MAGIC imputation is failed.")
                if show_error: redislogger.error(task_id, f"MAGIC imputation is failed: {e}")
        else: 
            redislogger.warning(task_id, "'MAGIC_imputed' layer already exists.")


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
        adata, counts, csv_path = load_anndata_to_csv(input, output, layer, show_error)
        if adata is None:
            redislogger.warning(task_id, f"File format is not supported: {input}")
        elif 'SAVER_imputed' not in adata.layers.keys(): 
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
            except Exception as e:
                redislogger.error(task_id, "SAVER imputation is failed.")
                if show_error: redislogger.error(task_id, f"SAVER imputation is failed: {e}")
        else: 
            redislogger.warning(task_id, "'SAVER_imputed' layer already exists.")

    return {'status': 'Success'}

    
    

        
            

