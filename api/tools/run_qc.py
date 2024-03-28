import os
import subprocess
import sys

from tools.qc.scanpy_qc import run_scanpy_qc
from tools.qc.dropkick_qc import run_dropkick_qc
from tools.qc.seurat_qc import run_seurat_qc
from tools.qc.scrublet_calls import predict_scrublet
# sys.path.append('..')
from tools.formating.formating import *
from config.celery_utils import get_output, benchmarks_output_path
from utils.unzip import unzip_file_if_compressed
from fastapi import HTTPException, status
from utils.redislogger import *
from tools.reduction.reduction import run_dimension_reduction, run_clustering
from utils.mongodb import generate_process_id, pp_results_exists, create_pp_results


def run_qc(task_id, ds:dict, random_state=0):
    results = []
    pp_results = []
    process_ids = []
    input_path = unzip_file_if_compressed(ds['userID'], ds['input'])
    methods = ds['methods']
    output = ds['output']
    adata_path = change_file_extension(input_path, 'h5ad')
    assay = ds['assay']
    assay_names = []
    md5 = get_md5(input_path)
    benchmarks_data = False
    if input_path is None:
        return None

    pp_stage = "Raw"
    process = "QC"
    if ds['qc_params']['max_genes'] == 20000:
        ds['qc_params']['max_genes'] = None
    if ds['qc_params']['n_pcs'] == 0:
        ds['qc_params']['n_pcs'] = None
    if ds['assay'] is None:
        ds['assay'] = 'RNA'
    
    parameters = ds['qc_params']
    print(type(parameters))
    redislogger.info(task_id, f"Using QC Parameters: {parameters}")
    
    # Get the absolute path for the given output
    output = get_output(output, ds['userID'], task_id) # Tools

    if methods is None: # Benchmarks, because benchmarks does not have method paramter
        benchmarks_data = True
        if input_path.endswith('.h5Seurat') or input_path.endswith('.h5seurat') or input_path.endswith('.rds') or input_path.endswith(".Robj"):
            methods = ["Seurat"]
        else:
            methods = ["scanpy"]
        output = benchmarks_output_path(input_path, task_id)

    # Get the absolute path for the given input
    # input = get_input_path(input, ds['userID'])
    # Get the absolute path for the given output
    methods = [x.upper() for x in methods if isinstance(x,str)]

    if "SCANPY" in methods or "DROPKICK" in methods:
        adata = load_anndata(input_path)
        # Scanpy QC
        if "SCANPY" in methods:
            method='scanpy'
            process_id = generate_process_id(md5, process, method, parameters, assay)
            qc_results = pp_results_exists(process_id)

            if qc_results is not None:
                redislogger.info(task_id, "Found existing pre-process results in database, skip Quality Control.")
            else:
                # Run Scanpy QC
                try:
                    output_path = get_output_path(output, ds['dataset'], method='scanpy')
                    # Check if the user only wants to run dimension reduction or clustering, then skip QC
                    # if do_qc:
                    redislogger.info(task_id, "Start scanpy QC...")
                    scanpy_results = run_scanpy_qc(adata, task_id, min_genes=parameters['min_genes'], max_genes=parameters['max_genes'], min_cells=parameters['min_cells'], target_sum=parameters['target_sum'], n_top_genes=parameters['n_top_genes'], expected_doublet_rate=parameters['doublet_rate'], regress_cell_cycle=parameters['regress_cell_cycle'])

                    # If the user only wants to run clustering, then skip dminension reduction
                    # if do_dimension:
                    redislogger.info(task_id, "Computing PCA, neighborhood graph, tSNE, UMAP, and 3D UMAP")
                    scanpy_results, msg = run_dimension_reduction(scanpy_results, n_neighbors=parameters['n_neighbors'], n_pcs=parameters['n_pcs'], random_state=random_state)
                    if msg is not None: redislogger.warning(task_id, msg)

                    # If the user only wants to run dminension reduction, then skip clustering
                    # if do_clustering:
                    redislogger.info(task_id, "Clustering the neighborhood graph.")
                    scanpy_results = run_clustering(scanpy_results, resolution=parameters['resolution'], random_state=random_state)
                    
                    redislogger.info(task_id, "Retrieving metadata and embeddings from AnnData object.")
                    qc_results = get_metadata_from_anndata(scanpy_results, pp_stage, process_id, process, method, parameters, adata_path)
                    redislogger.info(task_id, "Saving AnnData object.")
                    
                    scanpy_results.write_h5ad(output_path, compression='gzip')
                    scanpy_results = None
                    create_pp_results(qc_results)  # Insert pre-process results to database
                except Exception as e:
                    detail = f"Error during scanpy QC: {str(e)}"
                    redislogger.error(task_id, detail)
                    raise HTTPException(
                        status_code = status.HTTP_500_INTERNAL_SERVER_ERROR,
                        detail = detail
                    )
                 
            pp_results.append(qc_results)
            process_ids.append(process_id)

        # Dropkick QC
        if "DROPKICK" in methods:
            method='Dropkick'
            process_id = generate_process_id(md5, process, method, parameters,assay)
            qc_results = pp_results_exists(process_id)

            if qc_results is not None:
                redislogger.info(task_id, "Found existing pre-process results in database, skip Quality Control.")
            else:
                try:
                    redislogger.info(task_id, "Start Dropkick QC...")
                    dropkick_results = run_dropkick_qc(adata, task_id, n_neighbors=parameters['n_neighbors'], n_pcs=parameters['n_pcs'], resolution=parameters['resolution'], random_state=random_state)
                    
                    redislogger.info(task_id, "Computing PCA, neighborhood graph, tSNE, UMAP, and 3D UMAP")
                    dropkick_results, msg = run_dimension_reduction(dropkick_results, n_neighbors=parameters['n_neighbors'], n_pcs=parameters['n_pcs'], random_state=random_state)
                    if msg is not None: redislogger.warning(task_id, msg)

                    # If the user only wants to run dminension reduction, then skip clustering
                    # if do_clustering:
                    redislogger.info(task_id, "Clustering the neighborhood graph.")
                    dropkick_results = run_clustering(dropkick_results, resolution=parameters['resolution'], random_state=random_state)

                    redislogger.info(task_id, "Retrieving metadata and embeddings from AnnData object.")
                    qc_results = get_metadata_from_anndata(dropkick_results, pp_stage, process_id, process, method, parameters, adata_path)
                    redislogger.info(task_id, "Saving AnnData object.")
                    # output_path = get_output_path(output, ds['dataset'], method='dropkick')
                    output_path = "/usr/src/app/storage/kbcfh/Dataset2/Results/dropkick.h5ad"
                    redislogger.info(task_id, "output path")
                    redislogger.info(task_id, output_path)
                    dropkick_results.write_h5ad(output_path, compression='gzip')
                    dropkick_results = None
                    create_pp_results(qc_results) # Insert pre-process results to database
                except Exception as e:
                    detail = f"Error during Dropkick QC: {str(e)}"
                    redislogger.error(task_id, detail)
                    raise HTTPException(
                        status_code = status.HTTP_500_INTERNAL_SERVER_ERROR,
                        detail = detail
                    )
                
            pp_results.append(qc_results)
            process_ids.append(process_id)
            
        adata = None

    # Seurat QC
    if "SEURAT" in methods:
        method='Seurat'
        process_id = generate_process_id(md5, process, method, parameters,assay)
        qc_results = pp_results_exists(process_id)

        if qc_results is not None:
            redislogger.info(task_id, "Found existing pre-process results in database, skip Quality Control.")
        else:
            try:
                redislogger.info(task_id, "Start Seurat QC...")
                output_path = get_output_path(output, ds['dataset'], method='Seurat', format='Seurat')
                default_assay, assay_names, output, adata_path, adata, ddl_assay_names= run_seurat_qc(input_path, task_id, output=output_path, assay=ds['assay'], min_genes=parameters['min_genes'], max_genes=parameters['max_genes'], min_UMI_count=parameters['min_cells'], max_UMI_count=0, percent_mt_max=5, percent_rb_min=0, resolution=parameters['resolution'], dims=parameters['n_neighbors'], n_pcs=parameters['n_pcs'], doublet_rate=parameters['doublet_rate'], regress_cell_cycle=parameters['regress_cell_cycle'])
                
                if ddl_assay_names:
                    results.append({
                        "inputfile": input_path,
                        "default_assay": default_assay,
                        "assay_names": assay_names,
                        "ddl_assay_names": ddl_assay_names
                    })
                    return results
                
                redislogger.info(task_id, "Retrieving metadata and embeddings from AnnData object.")
                qc_results = get_metadata_from_anndata(adata, pp_stage, process_id, process, method, parameters, adata_path, seurat_path=output)
                create_pp_results(qc_results)  # Insert pre-process results to database
                adata = None         
            except Exception as e:
                detail = f"Error during Seurat QC: {str(e)}"
                redislogger.error(task_id, detail)
                raise HTTPException(
                    status_code = status.HTTP_500_INTERNAL_SERVER_ERROR,
                    detail = detail
                )
        if assay_names is None:
            assay_names = []
        
        pp_results.append(qc_results)
        process_ids.append(process_id)  

    # Bioconductor QC
    if "BIOCONDUCTOR" in methods:
        method='Bioconductor'
        process_id = generate_process_id(md5, process, method, parameters,assay)
        
        qc_results = pp_results_exists(process_id)

        if qc_results is not None:
            redislogger.info(task_id, "Found existing pre-process results in database, skip Quality Control.")
        else:
            try:
                redislogger.info(task_id, "Start Bioconductor QC...")
                output_path = get_output_path(output, ds['dataset'], method='Bioconductor', format='SingleCellExperiment')
                adata_path = get_output_path(output, ds['dataset'], method='Bioconductor', format='AnnData')
                report_path = get_report_path(ds['dataset'], output_path, "Bioconductor")

                # Get the absolute path of the current file
                current_file = os.path.abspath(__file__)

                # Construct the relative path to the desired file
                relative_path = os.path.join(os.path.dirname(current_file), 'qc', 'bioconductor_qc.Rmd')

                # Get the absolute path of the desired file
                bioconductor_path = os.path.abspath(relative_path)
                
                # bioconductor_path = os.path.abspath("qc/bioconductor_qc.Rmd")
                s = subprocess.call([f"R -e \"rmarkdown::render('{bioconductor_path}', params=list(dataset='{ds['dataset']}', input_path='{input_path}', idtype='{ds['idtype']}', colour_by='{parameters['colour_by']}', shape_by_1='{parameters['shape_by_1']}', shape_by_2='{parameters['shape_by_2'] }', output='{output_path}', adata_path='{adata_path}', output_format='SingleCellExperiment'), output_file='{report_path}')\""], shell = True)
                redislogger.info(task_id, s)

                if os.path.exists(adata_path):
                    redislogger.info(task_id, "Adding 3D UMAP to AnnData object.")
                    adata = load_anndata(adata_path)
                    sc.pp.neighbors(adata, n_neighbors=parameters['n_neighbors'], n_pcs=parameters['n_pcs'], random_state=0)
                    adata_3D = sc.tl.umap(adata, random_state=0, 
                                    init_pos="spectral", n_components=3, 
                                    copy=True, maxiter=None)
                    adata.obsm["X_umap_3D"] = adata_3D.obsm["X_umap"]
                    adata.write_h5ad(adata_path, compression='gzip')
                    adata_3D = None
                else:
                    raise ValueError("AnnData file does not exist due to the failure of Bioconductor QC.")
                
                redislogger.info(task_id, "Retrieving metadata and embeddings from AnnData object.")
                qc_results = get_metadata_from_anndata(adata, pp_stage, process_id, process, method, parameters, adata_path, sce_path=output_path)
                adata = None
                create_pp_results(qc_results)  # Insert pre-process results to database            
            except Exception as e:
                detail = f"Error during Bioconductor QC: {str(e)}"
                redislogger.error(task_id, detail)
                raise HTTPException(
                    status_code = status.HTTP_500_INTERNAL_SERVER_ERROR,
                    detail = detail
                )
        pp_results.append(qc_results)
        process_ids.append(process_id)    

    results.append({
        "task_id": task_id, 
        "inputfile": input_path,
        "default_assay": assay,
        "assay_names": assay_names,
        "md5": md5,
        "process_id": process_ids,
        "pp_results": pp_results,
        "message": "Quality control completed successfully."
    }) 

    return results