
import os
import subprocess
import sys
sys.path.append('..')
# from tools.formating.formating import *
from tools.evaluation.monitor import *
from tools.evaluation.integration import integration_metrics


def liger_integration(input, label, batch_key, benchmarksId, datasetId, task_type, cluster_key="leiden", species="mouse"):
    adata_int = None
    # Start monitoring
    monitor = Monitor(1)
    sys_info = monitor.get_sys_info()
    results = {}
    current_file = os.path.abspath(__file__)
    # Construct the relative path to the desired file
    relative_path = os.path.join(os.path.dirname(current_file), 'integration.Rmd')
    # Get the absolute path of the desired file
    rmd_path = os.path.abspath(relative_path)
    output = os.path.join(os.path.dirname(inputs[0]), f'{method}_integration')
    adata_path = input[0].replace(".h5ad", f"{method}_integration.h5ad")
    report_path = os.path.join(output, f'{method}_integration_report.html')
    s = subprocess.call([f"R -e \"rmarkdown::render('{rmd_path}', params=list(unique_id='{benchmarksId}', datasets='{datasetId}', inputs='{input}', output_folder='{output}', adata_path='{adata_path}', methods='LIGER'), output_file='{report_path}')\""], shell = True)
    
    # Stop monitoring
    time_points, cpu_usage, mem_usage, gpu_mem_usage = monitor.stop()

    current_date_and_time = datetime.now()

    if os.path.exists(adata_path):
        adata_int = load_anndata(adata_path)

    metrics_dict = integration_metrics(adata, adata_int, batch_key=batch_key, label_key=label, cluster_key=cluster_key, organism=species)

    sys_usage = {
            "sys_info": sys_info,
            "benchmarksId": benchmarksId,
            "datasetId": datasetId,
            "task_type": task_type,
            "tool": "Liger",
            "time_points": time_points,
            "cpu_usage": cpu_usage,
            "mem_usage": mem_usage,
            "gpu_mem_usage": gpu_mem_usage,
            "created_on": current_date_and_time
            }

    results["Liger"] = {**sys_usage, **metrics_dict}
    adata_int = None

    return results