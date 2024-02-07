from tools.formating.formating import load_anndata
from tools.visualization.plot import plot_bar, plot_line
from tools.benchmarks.clustering.scanpy import scanpy_clustering
from tools.benchmarks.clustering.scvi import scvi_clustering
from tools.benchmarks.clustering.seurat import seurat_clustering

def clustering_task(adata_path, task_label):
    # Load AnnData
    adata = load_anndata(adata_path)

    #static array to define the metrics evaluated for the clustering methods
    metrics = ['ARI', 'Silhouette', 'NMI']
    
    # Call scanpy_clustering method
    asw_scanpy, nmi_scanpy, ari_scanpy, time_points_scanpy, cpu_usage_scanpy, mem_usage_scanpy, gpu_mem_usage_scanpy = scanpy_clustering(adata, task_label)
    scanpy_results = {
        "asw_score": asw_scanpy,
        "nmi_score": nmi_scanpy,
        "ari_score": ari_scanpy,
        "time_points": time_points_scanpy,
        "cpu_usage": cpu_usage_scanpy,
        "mem_usage": mem_usage_scanpy,
        "gpu_mem_usage": gpu_mem_usage_scanpy
    }

    # Call seurat_clustering method
    asw_seurat, nmi_seurat, ari_seurat, time_points_seurat, cpu_usage_seurat, mem_usage_seurat, gpu_mem_usage_seurat = seurat_clustering(adata_path, task_label)

    seurat_results = {
        "asw_score": asw_seurat,
        "nmi_score": nmi_seurat,
        "ari_score": ari_seurat,
        "time_points": time_points_seurat,
        "cpu_usage": cpu_usage_seurat,
        "mem_usage": mem_usage_seurat,
        "gpu_mem_usage": gpu_mem_usage_seurat
    }
    
    # Call scvi_clustering method
    asw_scvi, nmi_scvi, ari_scvi, time_points_scvi, cpu_usage_scvi, mem_usage_scvi, gpu_mem_usage_scvi = scvi_clustering(adata_path, task_label)

    scvi_results = {
        "asw_score": asw_scvi,
        "nmi_score": nmi_scvi,
        "ari_score": ari_scvi,
        "time_points": time_points_scvi,
        "cpu_usage": cpu_usage_scvi,
        "mem_usage": mem_usage_scvi,
        "gpu_mem_usage": gpu_mem_usage_scvi
    }

    # Format x and y for the plot_bar function
    y_values = {
        'Scanpy': [ scanpy_results['ari_score'], scanpy_results['asw_score'], scanpy_results['nmi_score']],
        'Seurat': [seurat_results['ari_score'], seurat_results['asw_score'], seurat_results['nmi_score']],
        'scvi': [scvi_results['ari_score'], scvi_results['asw_score'], scvi_results['nmi_score']],
    }

    # # Call the plot_bar function
    bar_plot = plot_bar(x=metrics, y=y_values, title='Benchmarks')

    # Format x and y for the plot_line function
    y_values = {
        'Scanpy_CPU': scanpy_results['cpu_usage'],
        'Scanpy_Memory': scanpy_results['mem_usage'],
        'Scanpy_GPU': scanpy_results['gpu_mem_usage'],
        'Seurat_CPU': seurat_results['cpu_usage'],
        'Seurat_Memory': seurat_results['mem_usage'],
        'Seurat_GPU': seurat_results['gpu_mem_usage'],
        'scvi_CPU': scvi_results['cpu_usage'],
        'scvi_Memory': scvi_results['mem_usage'],
        'scvi_GPU': scvi_results['gpu_mem_usage'],
    }
    x_timepoints = []
    if len(scanpy_results['time_points']) > len(seurat_results['time_points']):
        x_timepoints = scanpy_results['time_points']
    else:
        x_timepoints =  seurat_results['time_points']
    
    if len(scvi_results['time_points']) > len(x_timepoints):
        x_timepoints = scvi_results['time_points']

    # Call the plot_line function with an empty array for x
    line_plot = plot_line(x=x_timepoints, y=y_values)
    
    result = {
        "adata_path": adata_path,
        "metrics": metrics,
        "Scanpy": scanpy_results,
        "Seurat": seurat_results,
        "scvi": scvi_results,
        "bar_plot": bar_plot,
        "line_plot": line_plot
    }

    return result