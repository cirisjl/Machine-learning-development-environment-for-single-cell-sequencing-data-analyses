import os
import subprocess
import sys
sys.path.append('..')
from tools.formating.formating import *
from tools.evaluation.monitor import *
from tools.evaluation.annotation import annotation_metrics


def singler_annotation(adata, adata_path, label, benchmarksId, datasetId, task_type, SingleR_ref, ref_path=None, species="mouse"):

    # report_path = get_report_path(dataset, output, "SAVER")
    report_path = adata_path.replace(".h5ad", "_report.html")
    output_folder = os.path.dirname(adata_path)
    
    # Get the absolute path of the current file
    current_file = os.path.abspath(__file__)

    # Construct the relative path to the desired file
    relative_path = os.path.join(os.path.dirname(current_file), 'singleR.Rmd')

    # Get the absolute path of the desired file
    singler_path = os.path.abspath(relative_path)
    # Start monitoring
    monitor = Monitor(1)
    sys_info = monitor.get_sys_info()
    results = {}

    if label is not None and ref_path is not None:
        s = subprocess.call([f"R -e \"rmarkdown::render('{singler_path}', params=list(unique_id='{benchmarksId}', dataset='{datasetId}', input='{adata_path}', output_folder='{output_folder}', species='{species}', reference='{SingleR_ref}', user_ref='{ref_path}', user_label='{label}'), output_file='{report_path}')\""], shell = True)
    else:
        s = subprocess.call([f"R -e \"rmarkdown::render('{singler_path}', params=list(unique_id='{benchmarksId}', dataset='{datasetId}', input='{adata_path}', output_folder='{output_folder}', species='{species}', reference='{SingleR_ref}'), output_file='{report_path}')\""], shell = True)
    
    # Stop monitoring
    time_points, cpu_usage, mem_usage, gpu_mem_usage = monitor.stop()

    current_date_and_time = datetime.now()

    csv_main = output_folder + "/results_main.csv"
    csv_fine = output_folder + "/results_fine.csv"
    csv_user = output_folder + "/results_user.csv"

    if os.path.exists(csv_main):
        df_main = pd.read_csv(csv_main, index_col=0)
        accuracy, f1_macro, f1_micro, f1_weighted = annotation_metrics(adata.obs[label], df_main['labels'])
        results[f"SingleR_main: {SingleR_ref}"] = {
                "sys_info": sys_info,
                "benchmarksId": benchmarksId,
                "datasetId": datasetId,
                "task_type": task_type,
                "tool": f"SingleR_main: {SingleR_ref}",
                "Accuracy": accuracy,
                "F1_macro": f1_macro,
                "F1_micro": f1_micro,
                "F1_weighted": f1_weighted,
                "time_points": time_points,
                "cpu_usage": cpu_usage,
                "mem_usage": mem_usage,
                "gpu_mem_usage": gpu_mem_usage,
                "created_on": current_date_and_time
            }
    
    if os.path.exists(csv_fine):
        df_fine = pd.read_csv(csv_fine, index_col=0)
        accuracy, f1_macro, f1_micro, f1_weighted = annotation_metrics(adata.obs[label], df_fine['labels'])
        results[f"SingleR_fine: {SingleR_ref}"] = {
                "sys_info": sys_info,
                "benchmarksId": benchmarksId,
                "datasetId": datasetId,
                "task_type": task_type,
                "tool": f"SingleR_fine: {SingleR_ref}",
                "Accuracy": accuracy,
                "F1_macro": f1_macro,
                "F1_micro": f1_micro,
                "F1_weighted": f1_weighted,
                "time_points": time_points,
                "cpu_usage": cpu_usage,
                "mem_usage": mem_usage,
                "gpu_mem_usage": gpu_mem_usage,
                "created_on": current_date_and_time
            }

    if os.path.exists(csv_user):
        df_user = pd.read_csv(csv_user, index_col=0)
        accuracy, f1_macro, f1_micro, f1_weighted = annotation_metrics(adata.obs[label], df_user['labels'])
        results[f"SingleR"] = {
                "sys_info": sys_info,
                "benchmarksId": benchmarksId,
                "datasetId": datasetId,
                "task_type": task_type,
                "tool": f"SingleR",
                "Accuracy": accuracy,
                "F1_macro": f1_macro,
                "F1_micro": f1_micro,
                "F1_weighted": f1_weighted,
                "time_points": time_points,
                "cpu_usage": cpu_usage,
                "mem_usage": mem_usage,
                "gpu_mem_usage": gpu_mem_usage,
                "created_on": current_date_and_time
            }

    return results