from tools.formating.formating import load_anndata, get_md5, clean_anndata
from tools.visualization.plot import plot_bar, plot_line
from benchmarks.multimodal_methods.multivi import multivi_multimodal
from utils.mongodb import generate_process_id, create_bm_results, benchmark_result_exists
from utils.redislogger import *
from datetime import datetime
import os


def multimodal_task(mdata_path, benchmarksId, datasetId, job_id, task_type='Multimodal'):
    redislogger.info(job_id, "Start running benchmarks for Multimodal task.")
    multimodal_results = []
    y_values = {}
    y_values_ur = {}
    x_timepoints = []
    md5 = get_md5(mdata_path)
    current_date_and_time = datetime.now()
    sys_info = None

    #static array to define the metrics evaluated for the Multimodal methods
    metrics = ["Mean Squared Error", "kNN Area Under the Curve"]
    
    # multivi
    try:
        redislogger.info(job_id, "Running MultiVI for Multimodal task.")
        process_id = generate_process_id(md5, task_type, 'MultiVI')
        multivi_results = benchmark_result_exists(process_id)

        if multivi_results is not None:
            redislogger.info(job_id, "Found existing MultiVI Benchmarks results in database, skip multivi.")
        else:
            # Call multivi method
            multivi_results = multivi_multimodal(mdata_path, benchmarksId, datasetId, task_type)
            create_bm_results(process_id, multivi_results)
            multimodal_results.append({'MultiVI': multivi_results})
            redislogger.info(job_id, "MultiVI Multimodal is done.")
        if len(multivi_results) > 0:
            for key, result in multivi_results.items():
                sys_info = result['sys_info']
                y_values[key] = [result['Mean Squared Error'], result['kNN Area Under the Curve']]
                y_values_ur['multivi_CPU'] = result['cpu_usage']
                y_values_ur['multivi_Memory'] = result['mem_usage']
                y_values_ur['multivi_GPU'] = result['gpu_usage']
                y_values_ur['multivi_GPU_Memory'] = result['gpu_mem_usage']
                x_timepoints = result['time_points']
                redislogger.info(job_id, f"{key}: Mean Squared Error: {result['Mean Squared Error']}, kNN Area Under the Curve: {result['kNN Area Under the Curve']}")

    except Exception as e:
        # Handle exceptions as needed
        redislogger.error(job_id, f"multivi Multimodal is failed: {e}")

    
    redislogger.info(job_id, "Creating bar plot for evaluation.")
    # Call the plot_bar function
    benchmarks_plot = plot_bar(x=metrics, y=y_values, title='Benchmarks: Multimodal')

    redislogger.info(job_id, "Creating line plot for computing resourses utilization rate.")
    # Call the plot_line function with an empty array for x
    utilization_plot = plot_line(x=x_timepoints, y=y_values_ur, sysinfo=sys_info)

    mdata = None # Release memory
    
    results = {
        # "mdata_path": mdata_path,
        "benchmarksId": benchmarksId,
        "datasetId": datasetId,
        "task_type": task_type,
        "metrics": metrics,
        "methods": multimodal_results,
        # "sys_info": sys_info,
        "benchmarks_plot": benchmarks_plot,
        "utilization_plot": utilization_plot
    }

    return results