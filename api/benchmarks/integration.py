from tools.formating.formating import load_anndata, get_md5, clean_anndata, get_scvi_path
from tools.visualization.plot import plot_bar, plot_line
from benchmarks.batch_integration_methods_methods.harmony import harmony_integration
from benchmarks.batch_integration_methods_methods.scvi import scvi_integration
from benchmarks.batch_integration_methods_methods.seurat import seurat_integration
from benchmarks.batch_integration_methods_methods.liger import liger_integration
from utils.mongodb import generate_process_id, create_bm_results, benchmark_result_exists
from utils.redislogger import *
from datetime import datetime
import os


def integration_task(adata_path, label, batch_key, benchmarksId, datasetId, job_id, species="mouse", task_type='Batch Integration'):
    redislogger.info(job_id, "Start running benchmarks for Batch Integration task.")
    integration_results = []
    y_values = {}
    y_values_ur = {}
    x_timepoints = []
    md5 = get_md5(adata_path)
    # Load AnnData
    adata = load_anndata(adata_path)
    scvi_path = get_scvi_path(adata_path)
    adata = clean_anndata(adata) # Remove outliers
    input_folder = os.path.dirname(adata_path)

    # Create input files for Seurat and Liger
    input = []
    for sample in adata.obs[batch_key].unique():
        adata_sub = adata[adata.obs[batch_key]==sample, :]
        adata_sub_path = f"input_folder/{sample}.h5ad"
        adata_sub.write_h5ad(adata_sub_path, compression='gzip')
        input.append(adata_sub_path)

    current_date_and_time = datetime.now()
    sys_info = None

    #static array to define the metrics evaluated for the Integration methods
    metrics = ['graph_conn', 'kBET', 'PCR_batch', 'ASW_label/batch', 'Biological Conservation']
    
    # Harmony
    try:
        redislogger.info(job_id, "Running Harmony for Batch Integration task.")
        process_id = generate_process_id(md5, task_type, 'Harmony', label)
        harmony_results = benchmark_result_exists(process_id)

        if harmony_results is not None:
            redislogger.info(job_id, "Found existing Harmony Benchmarks results in database, skip Harmony.")
        else:
            # Call Harmony method
            harmony_results = harmony_integration(adata, label, batch_key, benchmarksId, datasetId, task_type, species=species)
            create_bm_results(process_id, harmony_results)
            integration_results.append({'Harmony': harmony_results})
            redislogger.info(job_id, "Harmony Integration is done.")
        if len(harmony_results) > 0:
            for key, result in harmony_results.item
                sys_info = result['sys_info']
                y_values[key] = [result['graph_conn'], result['kBET'], result['PCR_batch'], result['ASW_label/batch'], result['Biological Conservation']]
                y_values_ur['Harmony_CPU'] = result['cpu_usage']
                y_values_ur['Harmony_Memory'] = result['mem_usage']
                y_values_ur['Harmony_GPU'] = result['gpu_mem_usage']
                x_timepoints = result['time_points']
                redislogger.info(job_id, f"{key}: Graph Connectivity: {result['graph_conn']}, kBET: {result['kBET']}, Principal component regression score: {result['PCR_batch']}, Batch ASW: {result['ASW_label/batch']}, Biological Conservation: {result['Biological Conservation']}")

    except Exception as e:
        # Handle exceptions as needed
        redislogger.error(job_id, f"Harmony Integration is failed: {e}")

    # scVI
    try:
        redislogger.info(job_id, "Running scVI for Batch Integration task.")
        process_id = generate_process_id(md5, task_type, 'scVI', label)
        scvi_results = benchmark_result_exists(process_id)

        if scvi_results is not None:
            redislogger.info(job_id, "Found existing scVI Benchmarks results in database, skip scVI.")
        else:
            # Call scVI method
            scvi_results = scvi_integration(adata, adata_path, label, batch_key, benchmarksId, datasetId, task_type, species=species)
            create_bm_results(process_id, scvi_results)
            integration_results.append({'scVI': scvi_results})
            redislogger.info(job_id, "scVI Integration is done.")
        if len(scvi_results) > 0:
            for key, result in scvi_results.item
                sys_info = result['sys_info']
                y_values[key] = [result['graph_conn'], result['kBET'], result['PCR_batch'], result['ASW_label/batch'], result['Biological Conservation']]
                y_values_ur['scVI_CPU'] = result['cpu_usage']
                y_values_ur['scVI_Memory'] = result['mem_usage']
                y_values_ur['scVI_GPU'] = result['gpu_mem_usage']
                x_timepoints = result['time_points']
                redislogger.info(job_id, f"{key}: Graph Connectivity: {result['graph_conn']}, kBET: {result['kBET']}, Principal component regression score: {result['PCR_batch']}, Batch ASW: {result['ASW_label/batch']}, Biological Conservation: {result['Biological Conservation']}")

    except Exception as e:
        # Handle exceptions as needed
        redislogger.error(job_id, f"scVI Integration is failed: {e}")
    
    # Seurat
    try:
        redislogger.info(job_id, "Running Seurat for Batch Integration task.")
        process_id = generate_process_id(md5, task_type, 'Seurat', label)
        seurat_results = benchmark_result_exists(process_id)

        if seurat_results is not None:
            redislogger.info(job_id, "Found existing Seurat Benchmarks results in database, skip Seurat.")
        else:
            # Call Seurat method
            seurat_results = seurat_integration(adata, label, batch_key, benchmarksId, datasetId, task_type, species=species)
            create_bm_results(process_id, seurat_results)
            integration_results.append({'Seurat': seurat_results})
            redislogger.info(job_id, "Seurat Integration is done.")
        if len(seurat_results) > 0:
            for key, result in seurat_results.item
                sys_info = result['sys_info']
                y_values[key] = [result['graph_conn'], result['kBET'], result['PCR_batch'], result['ASW_label/batch'], result['Biological Conservation']]
                y_values_ur['Seurat_CPU'] = result['cpu_usage']
                y_values_ur['Seurat_Memory'] = result['mem_usage']
                y_values_ur['Seurat_GPU'] = result['gpu_mem_usage']
                x_timepoints = result['time_points']
                redislogger.info(job_id, f"{key}: Graph Connectivity: {result['graph_conn']}, kBET: {result['kBET']}, Principal component regression score: {result['PCR_batch']}, Batch ASW: {result['ASW_label/batch']}, Biological Conservation: {result['Biological Conservation']}")

    except Exception as e:
        # Handle exceptions as needed
        redislogger.error(job_id, f"Seurat Integration is failed: {e}")

    # Liger
    try:
        redislogger.info(job_id, "Running Liger for Batch Integration task.")
        process_id = generate_process_id(md5, task_type, 'Liger', label)
        liger_results = benchmark_result_exists(process_id)

        if liger_results is not None:
            redislogger.info(job_id, "Found existing Liger Benchmarks results in database, skip Liger.")
        else:
            # Call Liger method
            liger_results = liger_integration(adata, label, batch_key, benchmarksId, datasetId, task_type, species=species)
            create_bm_results(process_id, liger_results)
            integration_results.append({'Liger': liger_results})
            redislogger.info(job_id, "Liger Integration is done.")
        if len(liger_results) > 0:
            for key, result in liger_results.item
                sys_info = result['sys_info']
                y_values[key] = [result['graph_conn'], result['kBET'], result['PCR_batch'], result['ASW_label/batch'], result['Biological Conservation']]
                y_values_ur['Liger_CPU'] = result['cpu_usage']
                y_values_ur['Liger_Memory'] = result['mem_usage']
                y_values_ur['Liger_GPU'] = result['gpu_mem_usage']
                x_timepoints = result['time_points']
                redislogger.info(job_id, f"{key}: Graph Connectivity: {result['graph_conn']}, kBET: {result['kBET']}, Principal component regression score: {result['PCR_batch']}, Batch ASW: {result['ASW_label/batch']}, Biological Conservation: {result['Biological Conservation']}")

    except Exception as e:
        # Handle exceptions as needed
        redislogger.error(job_id, f"Liger Integration is failed: {e}")
    
    redislogger.info(job_id, "Creating bar plot for evaluation.")
    # Call the plot_bar function
    benchmarks_plot = plot_bar(x=metrics, y=y_values, title='Benchmarks: Batch Integration')

    redislogger.info(job_id, "Creating line plot for computing resourses utilization rate.")
    # Call the plot_line function with an empty array for x
    utilization_plot = plot_line(x=x_timepoints, y=y_values_ur, sysinfo=sys_info)

    adata = None # Release memory
    
    results = {
        # "adata_path": adata_path,
        "benchmarksId": benchmarksId,
        "datasetId": datasetId,
        "task_type": task_type,
        "metrics": metrics,
        "methods": integration_results,
        # "sys_info": sys_info,
        "benchmarks_plot": benchmarks_plot,
        "utilization_plot": utilization_plot
    }

    return results