
import os
import subprocess
import sys
sys.path.append('..')
# from tools.formating.formating import *
from tools.evaluation.monitor import *
from tools.evaluation.integration import integration_metrics


def seurat_integration(input, label, batch_key, benchmarksId, datasetId, task_type, cluster_key="leiden", species="mouse"):
    # Start monitoring
    monitor = Monitor(1)
    sys_info = monitor.get_sys_info()
    results = {}
    current_file = os.path.abspath(__file__)
    # Construct the relative path to the desired file
    relative_path = os.path.join(os.path.dirname(current_file), 'integration.Rmd')
    # Get the absolute path of the desired file
    rmd_path = os.path.abspath(relative_path)
    s = subprocess.call([f"R -e \"rmarkdown::render('{rmd_path}', params=list(unique_id='{benchmarksId}', datasets='{datasetId}', inputs='{input}', output_folder='{output}', adata_path='{adata_path}', methods='SEURAT', dims={dims}, npcs={npcs}, resolution={resolution}, default_assay='{default_assay}'), output_file='{report_path}')\""], shell = True)
    
    # Stop monitoring
    time_points, cpu_usage, mem_usage, gpu_mem_usage = monitor.stop()

    current_date_and_time = datetime.now()

    metrics_dict = integration_metrics(adata, adata_int, batch_key=batch_key, label_key=label, cluster_key=cluster_key, organism=species)

    sys_usage = {
            "sys_info": sys_info,
            "benchmarksId": benchmarksId,
            "datasetId": datasetId,
            "task_type": task_type,
            "tool": "Seurat",
            "time_points": time_points,
            "cpu_usage": cpu_usage,
            "mem_usage": mem_usage,
            "gpu_mem_usage": gpu_mem_usage,
            "created_on": current_date_and_time
            }

    results["Seurat"] = {**sys_usage, **metrics_dict}

    return results