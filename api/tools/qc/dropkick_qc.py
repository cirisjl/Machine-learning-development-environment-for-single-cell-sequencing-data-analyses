import numpy as np
import pandas as pd
import scanpy as sc
import sklearn
import dropkick as dk
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

def dropkick_qc(adata):
    adata = dk.recipe_dropkick(adata, n_hvgs=None, X_final="raw_counts")
    
    # Run dropkick pipeline function
    qc_plt = dk.qc_summary(adata)

    adata_model = dk.dropkick(adata, n_jobs=5)

    score_plt = dk.score_plot(adata)

    dk.coef_inventory(adata)

    coef_plt = dk.coef_plot(adata)

    # Make a copy of our AnnData object, keeping all barcodes kept by dropkick, CellRanger, or EmptyDrops
    adata_filtered = adata[
        (adata.obs.dropkick_label=="True"),
        :].copy()
    
    # Here, we want to end up working with normalized, arcsinh-transformed counts 
    #   where genes are scaled to unit variance and zero-centered
    # we also set filter=True to remove any genes with zero total counts
    # we perform a variable gene selection for 2000 HVGs before further processing
    adata_filtered= dk.recipe_dropkick(adata_filtered, X_final="arcsinh_norm", filter=True, n_hvgs=2000, verbose=True)
    sc.pp.neighbors(adata_filtered, n_neighbors=30, random_state=1, n_pcs=10)
    sc.tl.leiden(adata_filtered)
    sc.tl.umap(adata_filtered, random_state=1)

    # plot results
    sc.pl.umap(
        adata_filtered,
        color=[
            "arcsinh_total_counts",
            "pct_counts_mito",
            "leiden",
            "dropkick_label",
            "dropkick_score",
        ],
        legend_fontsize="large",
        ncols=4,
    )

    # Save AnnData object
    # adata_filtered.write_h5ad(output, compression = 'gzip')
    # print("AnnData object for Dropkick QC is saved successfully")

    return adata_filtered
