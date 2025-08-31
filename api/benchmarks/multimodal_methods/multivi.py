import sys
sys.path.append('..')
# from tools.formating.formating import *
from tools.multimodal.MultiVI import run_multivi
from tools.evaluation.monitor import *
from tools.evaluation.multimodal import multimodal_metrics
from tools.formating.formating import *
from utils.redislogger import *
from datetime import datetime


def multivi_multimodal(mdata, benchmarksId, datasetId, task_type):
    # Start monitoring
    monitor = Monitor(1)
    sys_info = monitor.get_sys_info()
    results = {}
    mdata = run_multivi(mdata_path, rna_subset="rna_subset", atac_subset="atac_subset")
    
    # Stop monitoring
    time_points, cpu_usage, mem_usage, gpu_usage, gpu_mem_usage = monitor.stop()

    current_date_and_time = datetime.now()

    mse, area_under_curve = multimodal_metrics(mdata, aligned=mdata[mdata.obs["modality"]=='expression'].obsm['X_multivi'], mode2_aligned=mdata[mdata.obs["modality"]=='accessibility'].obsm['X_multivi'])

    results["MultiVI"] = {
                "sys_info": sys_info,
                "benchmarksId": benchmarksId,
                "datasetId": datasetId,
                "task_type": task_type,
                "tool": "MultiVI",
                "Mean Squared Error": mse,
                "kNN Area Under the Curve": area_under_curve,
                "time_points": time_points,
                "cpu_usage": cpu_usage,
                "mem_usage": mem_usage,
                "gpu_usage": gpu_usage,
                "gpu_mem_usage": gpu_mem_usage,
                "created_on": current_date_and_time
            }

    mdata = None

    return results