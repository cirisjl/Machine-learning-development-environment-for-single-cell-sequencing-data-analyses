from tools.formating.formating import load_anndata
from tools.visualization.plot import plot_bar, plot_line
from benchmarks.clustering_methods.scanpy import scanpy_clustering
from benchmarks.clustering_methods.scvi import scvi_clustering
from benchmarks.clustering_methods.seurat import seurat_clustering
from utils.redislogger import *

def clustering_task(adata_path, label, datasetId, task_id, task_type='clustering'):
    redislogger.info(task_id, "Start running benchmarks for clustering task.")
    clustering_results = []
    y_values = {}
    y_values_ur = {}
    x_timepoints = None
    # Load AnnData
    adata = load_anndata(adata_path)

    #static array to define the metrics evaluated for the clustering methods
    metrics = ['ARI', 'Silhouette', 'NMI']
    print(adata)

    print("run benchmarks")
    print(adata_path)
    print(label)
    print(task_type)
    print(datasetId)
    print(task_id)
    
    try:
        redislogger.info(task_id, "Running scanpy for clustering task.")
        # Call scanpy_clustering method
        asw_scanpy, nmi_scanpy, ari_scanpy, time_points_scanpy, cpu_usage_scanpy, mem_usage_scanpy, gpu_mem_usage_scanpy = scanpy_clustering(adata, label)
        scanpy_results = {
            "asw_score": asw_scanpy,
            "nmi_score": nmi_scanpy,
            "ari_score": ari_scanpy,
            "time_points": time_points_scanpy,
            "cpu_usage": cpu_usage_scanpy,
            "mem_usage": mem_usage_scanpy,
            "gpu_mem_usage": gpu_mem_usage_scanpy
        }
        y_values['scanpy'] = [scanpy_results['ari_score'], scanpy_results['asw_score'], scanpy_results['nmi_score']]
        y_values_ur['Scanpy_CPU'] = cpu_usage_scanpy
        y_values_ur['Scanpy_Memory'] = mem_usage_scanpy
        y_values_ur['Scanpy_GPU'] = gpu_mem_usage_scanpy
        x_timepoints = time_points_scanpy
        clustering_results.append({'scanpy': scanpy_results})

    except Exception as e:
        # Handle exceptions as needed
        redislogger.error(task_id, f"scanpy clustering is failed: {e}")

    try: 
        redislogger.info(task_id, "Running Seurat for clustering task.")
        # Call seurat_clustering method
        asw_seurat, nmi_seurat, ari_seurat, time_points_seurat, cpu_usage_seurat, mem_usage_seurat, gpu_mem_usage_seurat = seurat_clustering(adata_path, label)

        seurat_results = {
            "asw_score": asw_seurat,
            "nmi_score": nmi_seurat,
            "ari_score": ari_seurat,
            "time_points": time_points_seurat,
            "cpu_usage": cpu_usage_seurat,
            "mem_usage": mem_usage_seurat,
            "gpu_mem_usage": gpu_mem_usage_seurat
        }
        y_values['Seurat'] = [seurat_results['ari_score'], seurat_results['asw_score'], seurat_results['nmi_score']]
        y_values_ur['Seurat_CPU'] = cpu_usage_seurat
        y_values_ur['Seurat_Memory'] = mem_usage_seurat
        y_values_ur['Seurat_GPU'] = gpu_mem_usage_seurat
        if len(x_timepoints) < len(time_points_seurat):
            x_timepoints = time_points_seurat
        clustering_results.append({'Seurat': scanpy_results})

    except Exception as e:
        # Handle exceptions as needed
        redislogger.error(task_id, f"Seurat clustering is failed: {e}")
    
    try: 
        redislogger.info(task_id, "Running scvi for clustering task.")
        # Call scvi_clustering method
        asw_scvi, nmi_scvi, ari_scvi, time_points_scvi, cpu_usage_scvi, mem_usage_scvi, gpu_mem_usage_scvi = scvi_clustering(adata_path, label)

        scvi_results = {
            "asw_score": asw_scvi,
            "nmi_score": nmi_scvi,
            "ari_score": ari_scvi,
            "time_points": time_points_scvi,
            "cpu_usage": cpu_usage_scvi,
            "mem_usage": mem_usage_scvi,
            "gpu_mem_usage": gpu_mem_usage_scvi
        }
        y_values['scvi'] = [scvi_results['ari_score'], scvi_results['asw_score'], scvi_results['nmi_score']]
        y_values_ur['scvi_CPU'] = cpu_usage_scvi
        y_values_ur['scvi_Memory'] = mem_usage_scvi
        y_values_ur['scvi_GPU'] = gpu_mem_usage_scvi
        if len(x_timepoints) < len(time_points_scvi):
            x_timepoints = time_points_scvi
        clustering_results.append({'scvi': scvi_results})

    except Exception as e:
        # Handle exceptions as needed
        redislogger.error(task_id, f"scvi clustering is failed: {e}")
    
    redislogger.info(task_id, "Creating bar plot for evaluation.")
    # Call the plot_bar function
    bar_plot = plot_bar(x=metrics, y=y_values, title='Benchmarks')

    redislogger.info(task_id, "Creating line plot for computing resourses utilization rate.")
    # Call the plot_line function with an empty array for x
    line_plot = plot_line(x=x_timepoints, y=y_values_ur)
    
    results = {
        "adata_path": adata_path,
        "datasetId": datasetId,
        "task_type": task_type,
        "metrics": metrics,
        "results": clustering_results,
        "bar_plot": bar_plot,
        "line_plot": line_plot  
    }

    return results