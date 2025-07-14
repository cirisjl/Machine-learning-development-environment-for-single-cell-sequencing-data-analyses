from tools.formating.formating import load_anndata, get_md5, clean_anndata, get_scvi_path
from tools.visualization.plot import plot_bar, plot_line
from benchmarks.annotation_methods.celltypist import celltypist_annotation
from benchmarks.annotation_methods.scvi import scvi_annotation
from benchmarks.annotation_methods.singler import singler_annotation
from utils.mongodb import generate_process_id, create_bm_results, benchmark_result_exists
from utils.redislogger import *
from datetime import datetime
import os


def annotation_task(adata_path, label, benchmarksId, datasetId, job_id, celltypist_model=None, SingleR_ref=None, species="mouse", task_type='Cell Type Annotation'):
    redislogger.info(job_id, "Start running benchmarks for Cell Type Annotation task.")
    annotation_results = []
    y_values = {}
    y_values_ur = {}
    x_timepoints = []
    md5 = get_md5(adata_path)
    # Load AnnData
    adata = load_anndata(adata_path)
    scvi_path = get_scvi_path(adata_path)
    adata = clean_anndata(adata) # Remove outliers
    train_adata = adata[adata.obs.split_idx.str.contains('train'), :]
    ref_path = adata_path.replace(".h5ad", "_ref.h5ad")
    train_adata.write_h5ad(ref_path, compression='gzip')
    test_adata = adata[adata.obs.split_idx.str.contains('test'), :]
    current_date_and_time = datetime.now()
    sys_info = None

    #static array to define the metrics evaluated for the annotation methods
    metrics = ['Accuracy', 'F1_macro', 'F1_micro', 'F1_weighted']
    
    # CellTypist
    try:
        redislogger.info(job_id, "Running CellTypist for Cell Type Annotation task.")
        process_id = generate_process_id(md5, task_type, 'CellTypist', label)
        celltypist_results = benchmark_result_exists(process_id)

        if celltypist_results is not None:
            redislogger.info(job_id, "Found existing CellTypist Benchmarks results in database, skip CellTypist.")
        else:
            # Call CellTypist method
            celltypist_results = celltypist_annotation(adata, label, benchmarksId, datasetId, task_type, celltypist_model=celltypist_model, ref=train_adata, species=species)
            create_bm_results(process_id, celltypist_results)
            annotation_results.append({'CellTypist': celltypist_results})
            redislogger.info(job_id, "CellTypist annotation is done.")
        if len(celltypist_results) > 0:
            for key, result in celltypist_results.item
                sys_info = result['sys_info']
                y_values[key] = [result['Accuracy'], result['F1_macro'], result['F1_micro'], result['F1_weighted']]
                y_values_ur['CellTypist_CPU'] = result['cpu_usage']
                y_values_ur['CellTypist_Memory'] = result['mem_usage']
                y_values_ur['CellTypist_GPU'] = result['gpu_mem_usage']
                x_timepoints = result['time_points']
                redislogger.info(job_id, f"{key}: Accuracy: {result['Accuracy']}, F1_macro: {result['F1_macro']}, F1_micro: {result['F1_micro']}, F1_weighted: {result['F1_weighted']}")

    except Exception as e:
        # Handle exceptions as needed
        redislogger.error(job_id, f"CellTypist annotation is failed: {e}")

    # scVI
    try:
        redislogger.info(job_id, "Running scVI for Cell Type Annotation task.")
        process_id = generate_process_id(md5, task_type, 'scVI', label)
        scvi_results = benchmark_result_exists(process_id)

        if scvi_results is not None:
            redislogger.info(job_id, "Found existing scVI Benchmarks results in database, skip scVI.")
        else:
            # Call scVI method
            scvi_results = scvi_annotation(adata, label, benchmarksId, datasetId, task_type, ref=train_adata, species=species)
            create_bm_results(process_id, scvi_results)
            annotation_results.append({'scVI': scvi_results})
            redislogger.info(job_id, "scVI annotation is done.")
        if len(scvi_results) > 0:
            for key, result in scvi_results.item
                sys_info = result['sys_info']
                y_values[key] = [result['Accuracy'], result['F1_macro'], result['F1_micro'], result['F1_weighted']]
                y_values_ur['scVI_CPU'] = result['cpu_usage']
                y_values_ur['scVI_Memory'] = result['mem_usage']
                y_values_ur['scVI_GPU'] = result['gpu_mem_usage']
                x_timepoints = result['time_points']
                redislogger.info(job_id, f"{key}: Accuracy: {result['Accuracy']}, F1_macro: {result['F1_macro']}, F1_micro: {result['F1_micro']}, F1_weighted: {result['F1_weighted']}")

    except Exception as e:
        # Handle exceptions as needed
        redislogger.error(job_id, f"scVI annotation is failed: {e}")
    
    # SingleR
    try:
        redislogger.info(job_id, "Running SingleR for Cell Type Annotation task.")
        process_id = generate_process_id(md5, task_type, 'SingleR', label)
        singler_results = benchmark_result_exists(process_id)

        if singler_results is not None:
            redislogger.info(job_id, "Found existing SingleR Benchmarks results in database, skip SingleR.")
        else:
            # Call SingleR method
            singler_results = singler_annotation(adata, adata_path, label, benchmarksId, datasetId, task_type, SingleR_ref, ref_path=ref_path, species=species)
            create_bm_results(process_id, singler_results)
            annotation_results.append({'SingleR': singler_results})
            redislogger.info(job_id, "SingleR annotation is done.")
        if len(singler_results) > 0:
            for key, result in singler_results.item
                sys_info = result['sys_info']
                y_values[key] = [result['Accuracy'], result['F1_macro'], result['F1_micro'], result['F1_weighted']]
                y_values_ur['SingleR_CPU'] = result['cpu_usage']
                y_values_ur['SingleR_Memory'] = result['mem_usage']
                y_values_ur['SingleR_GPU'] = result['gpu_mem_usage']
                x_timepoints = result['time_points']
                redislogger.info(job_id, f"{key}: Accuracy: {result['Accuracy']}, F1_macro: {result['F1_macro']}, F1_micro: {result['F1_micro']}, F1_weighted: {result['F1_weighted']}")

    except Exception as e:
        # Handle exceptions as needed
        redislogger.error(job_id, f"SingleR annotation is failed: {e}")
    
    redislogger.info(job_id, "Creating bar plot for evaluation.")
    # Call the plot_bar function
    benchmarks_plot = plot_bar(x=metrics, y=y_values, title='Benchmarks: Cell Type Annotation')

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
        "methods": annotation_results,
        # "sys_info": sys_info,
        "benchmarks_plot": benchmarks_plot,
        "utilization_plot": utilization_plot
    }

    return results