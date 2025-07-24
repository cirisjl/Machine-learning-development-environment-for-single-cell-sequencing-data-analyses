from tools.formating.formating import load_anndata, load_anndata_to_csv, get_md5, clean_anndata, get_scvi_path
from tools.visualization.plot import plot_bar, plot_line
from benchmarks.imputation_methods.celltypist import celltypist_imputation
from benchmarks.imputation_methods.scvi import scvi_imputation
from utils.mongodb import generate_process_id, create_bm_results, benchmark_result_exists
from utils.redislogger import *
from datetime import datetime
import os


def imputation_task(adata_path, label, benchmarksId, datasetId, job_id, task_type='Imputation'):
    redislogger.info(job_id, "Start running benchmarks for Imputation task.")
    imputation_results = []
    y_values = {}
    y_values_ur = {}
    x_timepoints = []
    md5 = get_md5(adata_path)
    # Load AnnData
    csv_path = adata_path.replace(".h5ad", ".csv")
    adata, counts, csv_path = load_anndata_to_csv(adata_path, csv_path)
    current_date_and_time = datetime.now()
    sys_info = None

    #static array to define the metrics evaluated for the imputation methods
    metrics = ['MSE', 'Possion']
    
    # Magic
    try:
        redislogger.info(job_id, "Running Magic for Imputation task.")
        process_id = generate_process_id(md5, task_type, 'Magic', label)
        magic_results = benchmark_result_exists(process_id)

        if magic_results is not None:
            redislogger.info(job_id, "Found existing Magic Benchmarks results in database, skip Magic.")
        else:
            # Call Magic method
            magic_results = magic_imputation(adata, label, benchmarksId, datasetId, task_type, species=species)
            create_bm_results(process_id, magic_results)
            imputation_results.append({'Magic': magic_results})
            redislogger.info(job_id, "Magic imputation is done.")
        if len(magic_results) > 0:
            for key, result in magic_results.item
                sys_info = result['sys_info']
                y_values[key] = [result['MSE'], result['Possion']]
                y_values_ur['Magic_CPU'] = result['cpu_usage']
                y_values_ur['Magic_Memory'] = result['mem_usage']
                y_values_ur['Magic_GPU'] = result['gpu_mem_usage']
                x_timepoints = result['time_points']
                redislogger.info(job_id, f"{key}: MSE: {result['MSE']}, Possion: {result['Possion']}")

    except Exception as e:
        # Handle exceptions as needed
        redislogger.error(job_id, f"Magic imputation is failed: {e}")

   # Saver
    try:
        redislogger.info(job_id, "Running Saver for Imputation task.")
        process_id = generate_process_id(md5, task_type, 'Saver', label)
        saver_results = benchmark_result_exists(process_id)

        if saver_results is not None:
            redislogger.info(job_id, "Found existing Saver Benchmarks results in database, skip Saver.")
        else:
            # Call Saver method
            saver_results = saver_imputation(csv_path, label, benchmarksId, datasetId, task_type, species=species)
            create_bm_results(process_id, saver_results)
            imputation_results.append({'Saver': saver_results})
            redislogger.info(job_id, "Saver imputation is done.")
        if len(saver_results) > 0:
            for key, result in saver_results.item
                sys_info = result['sys_info']
                y_values[key] = [result['MSE'], result['Possion']]
                y_values_ur['Saver_CPU'] = result['cpu_usage']
                y_values_ur['Saver_Memory'] = result['mem_usage']
                y_values_ur['Saver_GPU'] = result['gpu_mem_usage']
                x_timepoints = result['time_points']
                redislogger.info(job_id, f"{key}: MSE: {result['MSE']}, Possion: {result['Possion']}")

    except Exception as e:
        # Handle exceptions as needed
        redislogger.error(job_id, f"Saver imputation is failed: {e}")
    
    redislogger.info(job_id, "Creating bar plot for evaluation.")
    # Call the plot_bar function
    benchmarks_plot = plot_bar(x=metrics, y=y_values, title='Benchmarks: Imputation')

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
        "methods": imputation_results,
        # "sys_info": sys_info,
        "benchmarks_plot": benchmarks_plot,
        "utilization_plot": utilization_plot
    }

    return results