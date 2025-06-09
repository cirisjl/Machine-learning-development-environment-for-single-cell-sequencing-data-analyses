import sys
import anndata as ad
from scib.metrics import metrics_all


def integration_metrics(adata, adata_int, ebed="X_pca", batch_key='batch', label_key='label_key', cluster_key = "leiden", organism="mouse"):
    """All metrics

    :Biological conservation:
        + HVG overlap :func:`~scib.metrics.hvg_overlap`
        + Cell type ASW :func:`~scib.metrics.silhouette`
        + Isolated label ASW :func:`~scib.metrics.isolated_labels`
        + Isolated label F1 :func:`~scib.metrics.isolated_labels`
        + NMI cluster/label :func:`~scib.metrics.nmi`
        + ARI cluster/label :func:`~scib.metrics.ari`
        + Cell cycle conservation :func:`~scib.metrics.cell_cycle`
        + cLISI (cell type Local Inverse Simpson's Index) :func:`~scib.metrics.clisi_graph`
        + Trajectory conservation :func:`~scib.metrics.trajectory_conservation`

    :Batch correction:
        + Graph connectivity :func:`~scib.metrics.graph_connectivity`
        + Batch ASW :func:`~scib.metrics.silhouette_batch`
        + Principal component regression :func:`~scib.metrics.pcr_comparison`
        + kBET (k-nearest neighbour batch effect test) :func:`~scib.metrics.kBET`
        + iLISI (integration Local Inverse Simpson's Index) :func:`~scib.metrics.ilisi_graph`

    :param adata: unintegrated, preprocessed anndata object
    :param adata_int: integrated anndata object
    :param batch_key: name of batch column in adata.obs and adata_int.obs
    :param label_key: name of biological label (cell type) column in adata.obs and adata_int.obs
    :param kwargs:
        Parameters to pass on to :func:`~scib.metrics.metrics` function:

            + ``embed``
            + ``cluster_key``
            + ``cluster_nmi``
            + ``nmi_method``
            + ``nmi_dir``
            + ``si_metric``
            + ``organism``
            + ``n_isolated``
            + ``subsample``
            + ``type_``
    """
    return metrics_all(adata, adata_int, ebed="X_pca", batch_key='batch', label_key='label_key', cluster_key = "leiden", organism="mouse")