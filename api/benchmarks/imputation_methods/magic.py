import sys
sys.path.append('..')
# from tools.formating.formating import *
from tools.imputation.MAGIC import magic_impute
from tools.evaluation.monitor import *
from tools.evaluation.imputation import imputation_metrics


def magic_imputation(adata, label, benchmarksId, datasetId, task_type, species="mouse"):
    # Start monitoring
    monitor = Monitor(1)
    sys_info = monitor.get_sys_info()
    results = {}
    counts = adata.X
    data_magic = magic_impute(counts, genes)
    adata.layers['MAGIC'] = data_magic
    
    # Stop monitoring
    time_points, cpu_usage, mem_usage, gpu_mem_usage = monitor.stop()

    current_date_and_time = datetime.now()

    mse, possion = imputation_metrics(adata, train="raw_counts", denoised=label, test="MAGIC")

    results["MAGIC"] = {
                "sys_info": sys_info,
                "benchmarksId": benchmarksId,
                "datasetId": datasetId,
                "task_type": task_type,
                "tool": "MAGIC",
                "MSE": mse,
                "Possion": possion,
                "time_points": time_points,
                "cpu_usage": cpu_usage,
                "mem_usage": mem_usage,
                "gpu_mem_usage": gpu_mem_usage,
                "created_on": current_date_and_time
            }

    adata = None

    return results