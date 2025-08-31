import sys
sys.path.append('..')
from tools.formating.formating import *
from tools.ccc.liana import run_liana_ccc
from tools.evaluation.monitor import *
from tools.evaluation.ccc import ccc_metrics
from datetime import datetime


def liana_ccc(adata, cell_type_label, benchmarksId, datasetId, task_type, species, ccc_pred='liana_res', ccc_target="ccc_target"):
    # Start monitoring
    monitor = Monitor(1)
    sys_info = monitor.get_sys_info()
    results = {}
    adata = run_liana_ccc(adata, cell_type_label=cell_type_label, species=pecies)
    
    # Stop monitoring
    time_points, cpu_usage, mem_usage, gpu_usage, gpu_mem_usage = monitor.stop()

    current_date_and_time = datetime.now()

    auc_score, oddsratio_score = ccc_metrics(adata, ccc_pred=ccc_pred, ccc_target=ccc_target)

    results["LIANA"] = {
                "sys_info": sys_info,
                "benchmarksId": benchmarksId,
                "datasetId": datasetId,
                "task_type": task_type,
                "tool": "LIANA",
                "Precision-recall AUC": auc_score,
                "Odds Ratio": oddsratio_score,
                "time_points": time_points,
                "cpu_usage": cpu_usage,
                "mem_usage": mem_usage,
                "gpu_usage": gpu_usage,
                "gpu_mem_usage": gpu_mem_usage,
                "created_on": current_date_and_time
            }

    adata = None

    return results