from tools.formating.formating import load_anndata, get_md5
from tools.visualization.plot import plot_bar, plot_line
from benchmarks.clustering_methods.scanpy import scanpy_clustering
from benchmarks.clustering_methods.scvi import scvi_clustering
from benchmarks.clustering_methods.seurat import seurat_clustering
from utils.mongodb import generate_process_id, create_bm_results, benchmark_result_exists
from utils.redislogger import *

def clustering_task(adata_path, label, datasetId, task_id, task_type='clustering'):
    redislogger.info(task_id, "Start running benchmarks for clustering task.")
    clustering_results = []
    y_values = {}
    y_values_ur = {}
    x_timepoints = None
    md5 = get_md5(adata_path)
    # Load AnnData
    adata = load_anndata(adata_path)

    #static array to define the metrics evaluated for the clustering methods
    metrics = ['ARI', 'Silhouette', 'NMI']
    
    # scanpy
    try:
        redislogger.info(task_id, "Running scanpy for clustering task.")
        process_id = generate_process_id(md5, task_type, 'scanpy', label)
        scanpy_results = benchmark_result_exists(process_id)

        if scanpy_results is not None:
            redislogger.info(task_id, "Found existing scanpy Benchmarks results in database, skip scanpy.")
        else:
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
            create_bm_results(process_id, scanpy_results)

        y_values['scanpy'] = [scanpy_results['ari_score'], scanpy_results['asw_score'], scanpy_results['nmi_score']]
        y_values_ur['Scanpy_CPU'] = scanpy_results['cpu_usage']
        y_values_ur['Scanpy_Memory'] = scanpy_results['mem_usage']
        y_values_ur['Scanpy_GPU'] = scanpy_results['gpu_mem_usage']
        x_timepoints = scanpy_results['time_points']
        clustering_results.append({'scanpy': scanpy_results})
        redislogger.info(task_id, "scanpy clustering is done.")
        redislogger.info(task_id, f"Silhouette: {scanpy_results['asw_score']}, NMI: {scanpy_results['nmi_score']}, ARI: {scanpy_results['ari_score']}")

    except Exception as e:
        # Handle exceptions as needed
        redislogger.error(task_id, f"scanpy clustering is failed: {e}")

    # Seurat
    try: 
        redislogger.info(task_id, "Running Seurat for clustering task.")
        process_id = generate_process_id(md5, task_type, 'Seurat', label)
        seurat_results = benchmark_result_exists(process_id)

        if seurat_results is not None:
            redislogger.info(task_id, "Found existing Seurat Benchmarks results in database, skip Seurat.")
        else:
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
            create_bm_results(process_id, seurat_results)

        y_values['Seurat'] = [seurat_results['ari_score'], seurat_results['asw_score'], seurat_results['nmi_score']]
        y_values_ur['Seurat_CPU'] = seurat_results['cpu_usage']
        y_values_ur['Seurat_Memory'] = seurat_results['mem_usage']
        y_values_ur['Seurat_GPU'] = seurat_results['gpu_mem_usage']
        if len(x_timepoints) < len(seurat_results['time_points']):
            x_timepoints = seurat_results['time_points']
        clustering_results.append({'Seurat': seurat_results})
        redislogger.info(task_id, "Seurat clustering is done.")
        redislogger.info(task_id, f"Silhouette: {seurat_results['asw_score']}, NMI: {seurat_results['nmi_score']}, ARI: {seurat_results['ari_score']}")

    except Exception as e:
        # Handle exceptions as needed
        redislogger.error(task_id, f"Seurat clustering is failed: {e}")
    
    # scVI
    try: 
        redislogger.info(task_id, "Running scvi for clustering task.")
        process_id = generate_process_id(md5, task_type, 'scVI', label)
        scvi_results = benchmark_result_exists(process_id)

        if scvi_results is not None:
            redislogger.info(task_id, "Found existing scVI Benchmarks results in database, skip scVI.")
        else:
            # Call scvi_clustering method
            asw_scvi, nmi_scvi, ari_scvi, time_points_scvi, cpu_usage_scvi, mem_usage_scvi, gpu_mem_usage_scvi = scvi_clustering(adata, label)

            scvi_results = {
                "asw_score": asw_scvi,
                "nmi_score": nmi_scvi,
                "ari_score": ari_scvi,
                "time_points": time_points_scvi,
                "cpu_usage": cpu_usage_scvi,
                "mem_usage": mem_usage_scvi,
                "gpu_mem_usage": gpu_mem_usage_scvi
            }
            create_bm_results(process_id, scvi_results)

        y_values['scvi'] = [scvi_results['ari_score'], scvi_results['asw_score'], scvi_results['nmi_score']]
        y_values_ur['scvi_CPU'] = scvi_results['cpu_usage']
        y_values_ur['scvi_Memory'] = scvi_results['mem_usage']
        y_values_ur['scvi_GPU'] = scvi_results['gpu_mem_usage']
        if len(x_timepoints) < len(scvi_results['time_points']):
            x_timepoints = scvi_results['time_points']
        clustering_results.append({'scvi': scvi_results})
        redislogger.info(task_id, "scvi clustering is done.")
        redislogger.info(task_id, f"Silhouette: {scvi_results['asw_score']}, NMI: {scvi_results['nmi_score']}, ARI: {scvi_results['ari_score']}")

    except Exception as e:
        # Handle exceptions as needed
        redislogger.error(task_id, f"scvi clustering is failed: {e}")
    
    redislogger.info(task_id, "Creating bar plot for evaluation.")
    # Call the plot_bar function
    bar_plot = plot_bar(x=metrics, y=y_values, title='Benchmarks')

    redislogger.info(task_id, "Creating line plot for computing resourses utilization rate.")
    # Call the plot_line function with an empty array for x
    line_plot = plot_line(x=x_timepoints, y=y_values_ur)

    adata = None # Release memory
    
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