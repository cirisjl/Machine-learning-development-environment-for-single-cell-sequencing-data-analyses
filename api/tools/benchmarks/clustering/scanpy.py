# Seed for reproducibility
import sys
sys.path.append('..')
# from tools.formating.formating import *
from tools.utils.utils import run_dimension_reduction
from tools.evaluation.monitor import *
from tools.evaluation.clustering import clustering_scores


def clustering(adata, labels, layer=None):
    # Start monitoring
    monitor = Monitor(1)

    adata = run_dimension_reduction(adata, layer=layer)
    
    # Stop monitoring
    time_points, cpu_usage, mem_usage, gpu_mem_usage = monitor.stop()

    if layer is None: layer = "X"
    # asw_score, nmi_score, ari_score = clustering_scores(adata.obs[labels], adata.obs["leiden"], adata.obsp['connectivities'])
    asw_score, nmi_score, ari_score = clustering_scores(adata.obs[labels], adata.obs["leiden"], adata.obsm[layer + '_umap'])

    return asw_score, nmi_score, ari_score, time_points, cpu_usage, mem_usage, gpu_mem_usage