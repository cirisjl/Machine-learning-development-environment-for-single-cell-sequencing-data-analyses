import numpy as np
import matplotlib.pyplot as pl
import scanpy as sc


def run_scanpy_trajectory(adata, color=["louvain"], n_neighbors=10, n_pcs=20, resolution=1.0):
    if adata is None:
        raise ValueError("Failed to load AnnData object.")
    
    if is_normalized(adata.X, 200) and not check_nonnegative_integers(adata.X):
        redislogger.info(unique_id, "adata.X is not raw counts.")
        if adata.raw.X is not None:
            redislogger.info(unique_id, "Use adata.raw.X instead. Copy adata.X to layer 'normalized_X'.")
            adata.layers["normalized_X"] = adata.X.copy()
            adata.X = adata.raw.X.copy()
        elif "raw_counts" in adata.layers.keys():
            redislogger.info(unique_id, "Use layer 'raw_counts' instead. Copy adata.X to layer 'normalized_X'.")
            adata.layers["normalized_X"] = adata.X.copy()
            adata.X = adata.layers['raw_counts'].copy()
        else:
            raise ValueError("scanpy trajectory only take raw counts, not normalized data.")

    adata.X = adata.X.astype("float64")

    sc.pp.recipe_zheng17(adata)
    sc.tl.pca(adata, svd_solver="arpack")
    # Denoising the graph
    sc.tl.diffmap(adata)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep="X_diffmap")
    sc.tl.louvain(adata, resolution=resolution)
    sc.tl.paga(adata, groups="louvain")
    sc.pl.paga(adata, color=color)
    # Recomputing the embedding using PAGA-initialization
    sc.tl.draw_graph(adata, init_pos="paga")

    sc.pl.paga_compare(
        adata,
        threshold=0.03,
        title="",
        right_margin=0.2,
        size=10,
        edge_width_scale=0.5,
        legend_fontsize=12,
        fontsize=12,
        frameon=False,
        edges=True,
        save=True,  # save figure to file figures/paga_compare.pdf
    )

    return adata
