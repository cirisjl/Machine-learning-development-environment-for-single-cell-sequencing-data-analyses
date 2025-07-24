import os
import subprocess
import sys
import shutil
from umap import UMAP
from tools.qc.scanpy_qc import run_scanpy_qc
# from tools.qc.dropkick_qc import run_dropkick_qc
from tools.qc.seurat_qc import run_seurat_qc
from tools.qc.scrublet_calls import predict_scrublet
# sys.path.append('..')
from tools.formating.formating import *
from config.celery_utils import get_output, benchmarks_output_path
from utils.unzip import unzip_file_if_compressed
from fastapi import HTTPException, status
from utils.redislogger import *
from tools.reduction.reduction import run_dimension_reduction, run_clustering
from utils.mongodb import generate_process_id, pp_result_exists, create_pp_results, upsert_jobs
from exceptions.custom_exceptions import CeleryTaskException
from datetime import datetime

def run_qc(job_id, ds:dict, fig_path=None, random_state=0):
    pp_results = []
    process_ids = []
    qc_output = []
    userID = ds['userID']
    input_path = unzip_file_if_compressed(job_id, ds['input'])
    datasetId = ds['datasetId']
    dataset = ds['dataset']
    output = ds['output']
    do_umap = ds['do_umap']
    do_cluster = ds['do_cluster']
    adata_path = change_file_extension(input_path, 'h5ad')
    assay_names = []
    md5 = get_md5(input_path)
    benchmarks_data = False
    description = "QC for Benchmarks"

    if datasetId is not None:
        description = f"QC for Benchmarks for {datasetId}"
    elif dataset is not None:
        description = f"QC for Benchmarks for {dataset}"

    if input_path is None:
        return None
    
    nCells = 0
    pp_stage = "Raw"
    process = "QC"
    parameters = ds['qc_params']
    assay = parameters['assay']
    if parameters['max_genes'] == 20000:
        parameters['max_genes'] = None
    if parameters['n_pcs'] == 0:
        parameters['n_pcs'] = None
    if assay is None:
        assay = 'RNA'
    skip_qc = parameters['skip_qc']
    
    methods = parameters['methods']
    method = None
    redislogger.info(job_id, f"Using QC Parameters: {parameters}")

    
    # Get the absolute path for the given output
    # output = get_output(output, ds['userID'], job_id) # Tools

    if methods is None and skip_qc == True: # Benchmarks, because benchmarks does not have method paramter
        redislogger.info(job_id, "Skip QC.")
        method = 'processed'
        adata = load_anndata(input_path)
        output = benchmarks_output_path(input_path)
        output_path = get_output_path(output, '', ds['dataset'])
        parameters.pop('min_genes')
        parameters.pop('max_genes')
        parameters.pop('min_cells')
        parameters.pop('target_sum')
        parameters.pop('doublet_rate')
        parameters.pop('regress_cell_cycle')
        parameters.pop('colour_by')
        parameters.pop('shape_by_1')
        parameters.pop('shape_by_2')

        process_id = generate_process_id(md5, process, method, parameters)
        qc_results = pp_result_exists(process_id)

        if qc_results is not None:
            redislogger.info(job_id, "Found existing pre-process results in database, skip Quality Control.")
            nCells = qc_results["nCells"]
            qc_output.append({'AnnData': qc_results["adata_path"]})
            adata_path = qc_results["adata_path"]
        else:
            if do_umap:
                redislogger.info(job_id, "Computing PCA, neighborhood graph, tSNE, UMAP, and 3D UMAP")
                adata, msg = run_dimension_reduction(adata, n_neighbors=parameters['n_neighbors'], n_pcs=parameters['n_pcs'], random_state=random_state, skip_if_exist=False)
                if msg is not None: redislogger.warning(job_id, msg)
            if do_cluster:
                redislogger.info(job_id, "Clustering the neighborhood graph.")
                adata = run_clustering(adata, resolution=parameters['resolution'], random_state=random_state, skip_if_exist=False)
            
            redislogger.info(job_id, "Retrieving metadata and embeddings from AnnData object.")
            qc_results = get_metadata_from_anndata(adata, pp_stage, process_id, process, method, parameters, md5, adata_path=output_path)
            nCells = qc_results["nCells"]
            redislogger.info(job_id, "Saving AnnData object.")
            adata.write_h5ad(output_path, compression='gzip')
            adata_path = output_path
            adata = None
            redislogger.info(job_id, qc_results['info'])
            qc_results['datasetId'] = datasetId
            qc_output.append({'AnnData': adata_path})
            create_pp_results(process_id, qc_results)  # Insert pre-process results to database
        process_ids.append(process_id)

    elif methods is None:
        benchmarks_data = True
        if input_path.endswith('.h5Seurat') or input_path.endswith('.h5seurat') or input_path.endswith('.rds') or input_path.endswith(".Robj"):
            methods = ["Seurat"]
        else:
            methods = ["scanpy"]
        
        upsert_jobs(
            {
                "job_id": job_id, 
                "description": description,
                "method": methods[0],
                "process": "Quality Control",
                "created_by": userID,
                "status": "Processing"
            }
        )
        output = benchmarks_output_path(input_path)
    else:
        upsert_jobs(
            {
                "job_id": job_id, 
                "created_by": userID,
                "status": "Processing"
            }
        )
    
    # Get the absolute path for the given input
    # input = get_input_path(input, ds['userID'])
    # Get the absolute path for the given output


    if skip_qc == False:
        methods = [x.upper() for x in methods if isinstance(x,str)]
        if ("SCANPY" in methods or "DROPKICK" in methods):
            # Scanpy QC
            if "SCANPY" in methods:
                method='scanpy'
                process_id = generate_process_id(md5, process, method, parameters)
                qc_results = pp_result_exists(process_id)

                if qc_results is not None:
                    redislogger.info(job_id, "Found existing pre-process results in database, skip Quality Control.")
                    nCells = qc_results["nCells"]
                    qc_output.append({'AnnData': qc_results["adata_path"]})
                    adata_path = qc_results["adata_path"]
                else:
                    output_path = None
                    if benchmarks_data:
                        output_path = get_output_path(output, '', ds['dataset'], method='scanpy')
                    else:
                        output_path = get_output_path(output, process_id, ds['dataset'], method='scanpy')
                    if os.path.exists(output_path): # If output exist from the last run, then just pick up it.
                        redislogger.info(job_id, "Output already exists, start from the result of the last run.")
                        scanpy_results = load_anndata(output_path)
                        if do_umap:
                            redislogger.info(job_id, "Computing PCA, neighborhood graph, tSNE, UMAP, and 3D UMAP")
                            scanpy_results, msg = run_dimension_reduction(scanpy_results, n_neighbors=parameters['n_neighbors'], n_pcs=parameters['n_pcs'], random_state=random_state)
                            if msg is not None: redislogger.warning(job_id, msg)                   
                        if do_cluster:
                            redislogger.info(job_id, "Clustering the neighborhood graph.")
                            scanpy_results = run_clustering(scanpy_results, resolution=parameters['resolution'], random_state=random_state)
                        
                        redislogger.info(job_id, "Retrieving metadata and embeddings from AnnData object.")
                        qc_results = get_metadata_from_anndata(scanpy_results, pp_stage, process_id, process, method, parameters, md5, adata_path=output_path, n_top_genes=parameters['n_top_genes'])
                        nCells = qc_results["nCells"]
                        redislogger.info(job_id, "Saving AnnData object.")
                        scanpy_results.X = scanpy_results.layers["raw_counts"].copy()
                        scanpy_results.write_h5ad(output_path, compression='gzip')
                        qc_output.append({'AnnData': output_path})
                        adata_path = output_path
                        scanpy_results = None
                        redislogger.info(job_id, qc_results['info'])
                        qc_results['datasetId'] = datasetId
                        create_pp_results(process_id, qc_results)  # Insert pre-process results to database
                    else:
                        # Run Scanpy QC 
                        try:
                            adata = load_anndata(input_path)
                            redislogger.info(job_id, "Start scanpy QC...")
                            scanpy_results = run_scanpy_qc(adata, job_id, min_genes=parameters['min_genes'], max_genes=parameters['max_genes'], min_cells=parameters['min_cells'], target_sum=parameters['target_sum'], n_top_genes=parameters['n_top_genes'], expected_doublet_rate=parameters['doublet_rate'], regress_cell_cycle=parameters['regress_cell_cycle'])
                            scanpy_results.write_h5ad(output_path, compression='gzip')

                            if do_umap:
                                redislogger.info(job_id, "Computing PCA, neighborhood graph, tSNE, UMAP, and 3D UMAP")
                                scanpy_results, msg = run_dimension_reduction(scanpy_results, n_neighbors=parameters['n_neighbors'], n_pcs=parameters['n_pcs'], random_state=random_state)
                                if msg is not None: redislogger.warning(job_id, msg)
                            if do_cluster:
                                redislogger.info(job_id, "Clustering the neighborhood graph.")
                                scanpy_results = run_clustering(scanpy_results, resolution=parameters['resolution'], random_state=random_state)
                            
                            redislogger.info(job_id, "Retrieving metadata and embeddings from AnnData object.")
                            qc_results = get_metadata_from_anndata(scanpy_results, pp_stage, process_id, process, method, parameters, md5, adata_path=output_path)
                            nCells = qc_results["nCells"]
                            redislogger.info(job_id, "Saving AnnData object.")
                            scanpy_results.X = scanpy_results.layers["raw_counts"].copy()
                            scanpy_results.write_h5ad(output_path, compression='gzip')
                            qc_output.append({'AnnData': output_path})
                            adata_path = output_path
                            scanpy_results = None
                            redislogger.info(job_id, qc_results['info'])
                            qc_results['datasetId'] = datasetId
                            
                            create_pp_results(process_id, qc_results)  # Insert pre-process results to database
                        except Exception as e:
                            detail = f"Error during scanpy QC: {str(e)}"
                            upsert_jobs(
                                {
                                    "job_id": job_id, 
                                    "results": detail,
                                    "completed_on": datetime.now(),
                                    "status": "Failure"
                                }
                            )
                            redislogger.error(job_id, detail)
                            raise CeleryTaskException(detail)
                    
                pp_results.append(qc_results)
                process_ids.append(process_id)

            # # Dropkick QC
            # if "DROPKICK" in methods:
            #     method='Dropkick'
            #     process_id = generate_process_id(md5, process, method, parameters)
            #     qc_results = pp_result_exists(process_id)

            #     if qc_results is not None:
            #         redislogger.info(job_id, "Found existing pre-process results in database, skip Quality Control.")
            #         nCells = qc_results["nCells"]
            #         qc_output.append({'AnnData': qc_results["adata_path"]})
            #         adata_path = qc_results["adata_path"]
            #     else:
            #         output_path = get_output_path(output, process_id, ds['dataset'], method='dropkick')
            #         if os.path.exists(output_path): # If output exist from the last run, then just pick up it.
            #             redislogger.info(job_id, "Output already exists, start from the result of the last run.")
            #             dropkick_results = load_anndata(output_path)
            
            #             if do_umap:
            #                 redislogger.info(job_id, "Computing PCA, neighborhood graph, tSNE, UMAP, and 3D UMAP")
            #                 dropkick_results.X = dropkick_results.layers["arcsinh_norm"].copy()
            #                 dropkick_results, msg = run_dimension_reduction(dropkick_results, n_neighbors=parameters['n_neighbors'], n_pcs=parameters['n_pcs'], random_state=random_state)
            #                 if msg is not None: redislogger.warning(job_id, msg)
            #             if do_cluster:
            #                 redislogger.info(job_id, "Clustering the neighborhood graph.")
            #                 dropkick_results = run_clustering(dropkick_results, resolution=parameters['resolution'], random_state=random_state)

            #             redislogger.info(job_id, "Retrieving metadata and embeddings from AnnData object.")
            #             qc_results = get_metadata_from_anndata(dropkick_results, pp_stage, process_id, process, method, parameters, md5, adata_path=output_path)
            #             nCells = qc_results["nCells"]
            #             redislogger.info(job_id, "Saving AnnData object.")
            #             dropkick_results.X = dropkick_results.layers["raw_counts"].copy()
            #             dropkick_results.write_h5ad(output_path, compression='gzip')
            #             qc_output.append({'AnnData': output_path})
            #             adata_path = output_path
            #             dropkick_results = None
            #             redislogger.info(job_id, qc_results['info'])
            #             qc_results['datasetId'] = datasetId
            #             create_pp_results(process_id, qc_results) # Insert pre-process results to database
            #         else:
            #             try:
            #                 redislogger.info(job_id, "Start Dropkick QC...")
            #                 adata = load_anndata(input_path)
            #                 dropkick_results = run_dropkick_qc(adata, job_id, n_neighbors=parameters['n_neighbors'], n_pcs=parameters['n_pcs'], resolution=parameters['resolution'], n_hvg=parameters['n_top_genes'], random_state=random_state)
            #                 dropkick_results.write_h5ad(output_path, compression='gzip')

            #                 if do_umap:
            #                     redislogger.info(job_id, "Computing PCA, neighborhood graph, tSNE, UMAP, and 3D UMAP")
            #                     dropkick_results.X = dropkick_results.layers["arcsinh_norm"].copy()
            #                     dropkick_results, msg = run_dimension_reduction(dropkick_results, n_neighbors=parameters['n_neighbors'], n_pcs=parameters['n_pcs'], random_state=random_state)
            #                     if msg is not None: redislogger.warning(job_id, msg)
            #                 if do_cluster:
            #                     redislogger.info(job_id, "Clustering the neighborhood graph.")
            #                     dropkick_results = run_clustering(dropkick_results, resolution=parameters['resolution'], random_state=random_state)

            #                 redislogger.info(job_id, "Retrieving metadata and embeddings from AnnData object.")
            #                 qc_results = get_metadata_from_anndata(dropkick_results, pp_stage, process_id, process, method, parameters, md5, adata_path=output_path)
            #                 nCells = qc_results["nCells"]
            #                 redislogger.info(job_id, "Saving AnnData object.")
            #                 # redislogger.info(job_id, "output path")
            #                 # redislogger.info(job_id, output_path)
            #                 dropkick_results.X = dropkick_results.layers["raw_counts"].copy()
            #                 dropkick_results.write_h5ad(output_path, compression='gzip')
            #                 adata_path = output_path
            #                 qc_output.append({'AnnData': output_path})
            #                 dropkick_results = None
            #                 redislogger.info(job_id, qc_results['info'])
            #                 qc_results['datasetId'] = datasetId
            #                 create_pp_results(process_id, qc_results) # Insert pre-process results to database
            #             except Exception as e:
            #                 detail = f"Error during Dropkick_QC: {str(e)}"
            #                 redislogger.error(job_id, detail)
            #                 upsert_jobs(
            #                     {
            #                         "job_id": job_id, 
            #                         "results": detail,
            #                         "completed_on": datetime.now(),
            #                         "status": "Failure"
            #                     }
            #                 )
            #                 os.remove(output_path)
            #                 raise CeleryTaskException(detail)
                    
            #     pp_results.append(qc_results)
            #     process_ids.append(process_id)
                
            # adata = None

        # Seurat QC
        if "SEURAT" in methods:
            method='Seurat'
            process_id = generate_process_id(md5, process, method, parameters)
            qc_results = pp_result_exists(process_id)

            if qc_results is not None:
                redislogger.info(job_id, "Found existing pre-process results in database, skip Quality Control.")
                nCells = qc_results["nCells"]
                qc_output.append({'AnnData': qc_results["adata_path"]})
                adata_path = qc_results["adata_path"]
            else:
                output_path = None
                if benchmarks_data:
                    output_path = get_output_path(output, '', ds['dataset'], method='Seurat', format='Seurat')
                else:
                    output_path = get_output_path(output, process_id, ds['dataset'], method='Seurat', format='Seurat')
                try:     
                    default_assay, assay_names, output_path, adata_path, adata, ddl_assay_names = run_seurat_qc(input_path, job_id, output=output_path, assay=assay, min_genes=parameters['min_genes'], max_genes=parameters['max_genes'], min_UMI_count=parameters['min_cells'], max_UMI_count=0, percent_mt_max=5, percent_rb_min=0, resolution=parameters['resolution'], dims=parameters['n_neighbors'], doublet_rate=parameters['doublet_rate'], n_hvg=parameters['n_top_genes'], regress_cell_cycle=parameters['regress_cell_cycle'])
                    
                    if ddl_assay_names:
                        results = {
                            "job_id": job_id, 
                            "inputfile": input_path,
                            "default_assay": default_assay,
                            "assay_names": assay_names,
                            "ddl_assay_names": ddl_assay_names
                        }
                        # upsert_jobs(results)
                        return results

                    if adata is not None:
                        if do_umap:
                            redislogger.info(job_id, "Computing PCA, neighborhood graph, tSNE, UMAP, and 3D UMAP")
                            adata.X = adata.layers["scale.data"].copy()
                            adata, msg = run_dimension_reduction(adata, n_neighbors=parameters['n_neighbors'], n_pcs=parameters['n_pcs'], random_state=random_state)
                            if msg is not None: redislogger.warning(job_id, msg)
                        if do_cluster:
                            redislogger.info(job_id, "Clustering the neighborhood graph.")
                            adata = run_clustering(adata, resolution=parameters['resolution'], random_state=random_state)
                        adata.X = adata.layers["raw_counts"].copy()

                        redislogger.info(job_id, "Retrieving metadata and embeddings from AnnData object.")
                        qc_results = get_metadata_from_anndata(adata, pp_stage, process_id, process, method, parameters, md5, adata_path=adata_path, seurat_path=output_path)
                        nCells = qc_results["nCells"]
                        redislogger.info(job_id, qc_results['info'])
                        adata.write_h5ad(adata_path, compression='gzip')
                        if os.path.exists(adata_path): qc_output.append({'Anndata': adata_path})
                        if os.path.exists(output_path): qc_output.append({'Seurat': output_path})
                        qc_results['datasetId'] = datasetId
                        create_pp_results(process_id, qc_results)  # Insert pre-process results to database
                        adata = None         
                except Exception as e:
                    detail = f"Error during Seurat QC: {str(e)}"
                    redislogger.error(job_id, detail)
                    upsert_jobs(
                        {
                            "job_id": job_id, 
                            "results": detail,
                            "completed_on": datetime.now(),
                            "status": "Failure"
                        }
                    )
                    os.remove(output_path)
                    os.remove(adata_path)
                    raise CeleryTaskException(detail)
            if assay_names is None:
                assay_names = []
            
            pp_results.append(qc_results)
            process_ids.append(process_id)  

        # Bioconductor QC
        if "BIOCONDUCTOR" in methods:
            method='Bioconductor'
            process_id = generate_process_id(md5, process, method, parameters)
            qc_results = pp_result_exists(process_id)

            if qc_results is not None:
                redislogger.info(job_id, "Found existing pre-process results in database, skip Quality Control.")
                nCells = qc_results["nCells"]
                qc_output.append({"AnnData": qc_results["adata_path"]})
                adata_path = qc_results["adata_path"]
            else:
                output_path = get_output_path(output, process_id, ds['dataset'], method='Bioconductor', format='SingleCellExperiment')
                adata_path = get_output_path(output, process_id, ds['dataset'], method='Bioconductor', format='AnnData')
                report_path = get_report_path(ds['dataset'], output_path, "Bioconductor")
                try:
                    redislogger.info(job_id, "Start Bioconductor QC...")

                    # Get the absolute path of the current file
                    current_file = os.path.abspath(__file__)

                    # Construct the relative path to the desired file
                    relative_path = os.path.join(os.path.dirname(current_file), 'qc', 'bioconductor_qc.Rmd')

                    # Get the absolute path of the desired file
                    bioconductor_path = os.path.abspath(relative_path)
                    
                    # bioconductor_path = os.path.abspath("qc/bioconductor_qc.Rmd")
                    s = subprocess.call([f"R -e \"rmarkdown::render('{bioconductor_path}', params=list(dataset='{ds['dataset']}', input_path='{input_path}', idtype='{ds['idtype']}', colour_by='{parameters['colour_by']}', shape_by_1='{parameters['shape_by_1']}', shape_by_2='{parameters['shape_by_2'] }', output='{output_path}', adata_path='{adata_path}', output_format='SingleCellExperiment', n_hvg={parameters['n_top_genes']}), output_file='{report_path}')\""], shell = True)
                    # redislogger.info(job_id, s)

                    if os.path.exists(adata_path):
                        redislogger.info(job_id, "Adding 3D UMAP to AnnData object.")
                        adata = load_anndata(adata_path)
                        if 'subset' in adata.var.columns:
                            adata.var = adata.var.rename(columns={'subset': 'highly_variable'})
                        adata.layers["raw_counts"] = adata.X.copy()
                        sc.pp.normalize_total(adata, target_sum=parameters['target_sum'])
                        sc.pp.log1p(adata)
                        if do_umap:
                            redislogger.info(job_id, "Computing PCA, neighborhood graph, tSNE, UMAP, and 3D UMAP")
                            adata, msg = run_dimension_reduction(adata, n_neighbors=parameters['n_neighbors'], n_pcs=parameters['n_pcs'], random_state=random_state)
                            if msg is not None: redislogger.warning(job_id, msg)
                        if do_cluster:
                            redislogger.info(job_id, "Clustering the neighborhood graph.")
                            adata = run_clustering(adata, resolution=parameters['resolution'], random_state=random_state)
                        adata.X = adata.layers["raw_counts"].copy()
                        # sc.pp.neighbors(adata, n_neighbors=parameters['n_neighbors'], n_pcs=parameters['n_pcs'], random_state=0)
                        # umap_3d = UMAP(n_components=3, init='random', random_state=0)
                        # adata.obsm["X_umap_3D"] = umap_3d.fit_transform(adata.obsm['X_pca'])
                        # adata_3D = sc.tl.umap(adata, random_state=0, 
                        #                 init_pos="spectral", n_components=3, 
                        #                 copy=True, maxiter=None)
                        # adata.obsm["X_umap_3D"] = adata_3D.obsm["X_umap"]
                        adata.write_h5ad(adata_path, compression='gzip')
                    else:
                        raise ValueError("AnnData file does not exist due to the failure of Bioconductor QC.")
                    
                    redislogger.info(job_id, "Retrieving metadata and embeddings from AnnData object.")
                    qc_results = get_metadata_from_anndata(adata, pp_stage, process_id, process, method, parameters, md5, adata_path=adata_path, sce_path=output_path)
                    nCells = qc_results["nCells"]
                    if os.path.exists(adata_path): qc_output.append({'Anndata': adata_path})
                    if os.path.exists(output_path): qc_output.append({'SingleCellExperiment': output_path})
                    if os.path.exists(report_path): qc_output.append({'Report': report_path})
                    
                    adata = None
                    redislogger.info(job_id, qc_results['info'])
                    qc_results['datasetId'] = datasetId
                    create_pp_results(process_id, qc_results)  # Insert pre-process results to database            
                except Exception as e:
                    detail = f"Error during Bioconductor QC: {str(e)}"
                    redislogger.error(job_id, detail)
                    upsert_jobs(
                        {
                            "job_id": job_id, 
                            "results": detail,
                            "completed_on": datetime.now(),
                            "status": "Failure"
                        }
                    )
                    os.remove(output_path)
                    os.remove(adata_path)
                    os.remove(report_path)
                    raise CeleryTaskException(detail)
                    
            pp_results.append(qc_results)
            process_ids.append(process_id)    

    results = {
        "nCells": nCells,
        "owner": userID,
        "method": method,
        "output": qc_output,
        "adata_path": adata_path,
        "default_assay": assay,
        "assay_names": assay_names,
        "md5": md5,
        "process_ids": process_ids
    }

    upsert_jobs(
        {
            "job_id": job_id, 
            "datasetId": datasetId,
            "process_ids": process_ids,
            "nCells": nCells,
            "output": qc_output,
            "default_assay": assay,
            "assay_names": assay_names,
            "results": results,
            "completed_on": datetime.now(),
            "status": "Success"
        }
    )

    # Remove oringal input files after QC for Benchmarks. Input files will be removed when submiting metadata form
    # if benchmarks_data:
    #     if os.path.isdir(input_path):
    #         shutil.rmtree(input_path) 
    #     else:
    #         os.remove(input_path)

    return results