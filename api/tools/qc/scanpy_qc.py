import numpy as np
import pandas as pd
import scanpy as sc
import scrublet as scr
import warnings
warnings.filterwarnings('ignore')
import sklearn
import scipy
# sys.path.append('..')
from scipy.stats import median_abs_deviation
from tools.formating.formating import is_normalized
from scipy.sparse import csr_matrix
sc.settings.verbosity=3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')


def run_scanpy_qc(adata, min_genes=200, min_cells=3, target_sum=1e4, regress_cell_cycle=False):
        if adata is None:
            raise ValueError("The input is None.")
        
        if is_normalized(adata.X, min_genes):
            if adata.raw.X is not None:
                adata.layers["normalized_X"] = adata.X
                adata.X = adata.raw.X
            elif "raw_counts" in adata.layers.keys():
                adata.layers["normalized_X"] = adata.X
                adata.X = adata.layers['raw_counts']

        if is_normalized(adata.X, min_genes):
            raise ValueError("Scanpy QC only take raw counts, not normalized data.")
        
        adata.var_names_make_unique()

        # Filtering low quality reads
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=min_cells)
        # mitochondrial genes
        adata.var['mt']=adata.var_names.str.startswith('MT-')
        # ribosomal genes
        adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
        # hemoglobin genes.
        adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]"))
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True)

        adata.obs["outlier"] = (
            is_outlier(adata, "log1p_total_counts", 5)
            | is_outlier(adata, "log1p_n_genes_by_counts", 5)
            | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
        )
        adata.obs.outlier.value_counts()

        adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (
            adata.obs["pct_counts_mt"] > 8
        )
        adata.obs.mt_outlier.value_counts()

        print(f"Total number of cells: {adata.n_obs}")
        adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()

        print(f"Number of cells after filtering of low quality cells: {adata.n_obs}")

        # adata=adata[adata.obs.n_genes_by_counts < 2500, :]
        # adata=adata[adata.obs.pct_counts_mt < 5, :]

        scrub = scr.Scrublet(adata.X, expected_doublet_rate = 0.076)
        adata.obs['doublet_scores'], adata.obs['predicted_doublets'] = scrub.scrub_doublets(min_counts=2, min_cells=3, 
                                                                min_gene_variability_pctl=85, n_prin_comps=30)
        adata.obs['predicted_doublets'].value_counts()
        # adata=adata[adata.obs.predicted_doublets=="False", :]

        if adata.raw is None:
            adata.raw = adata
        else: 
            adata.layers["raw_counts"] = adata.X
        
        sc.pp.normalize_total(adata, target_sum=target_sum)

        sc.pp.log1p(adata)

        sc.pp.highly_variable_genes(adata)

        adata=adata[:, adata.var.highly_variable] # Do the filtering
        sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
        sc.pp.scale(adata, max_value=10)

        # Regress both S score and G2M score for cell cycle
        if(regress_cell_cycle):
             adata = regress_cell_cycle(adata)

        if isinstance(adata.X, np.ndarray):
            adata.X = csr_matrix(adata.X)

        adata.layers["log10k"] = adata.X

        # return adata, output
        return adata


def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier


def regress_cell_cycle(adata):
    # Load cell cycle genes defined in Tirosh et al, 2015. It is a list of 97 genes, represented by their gene symbol.
    cell_cycle_genes = [x.strip() for x in open(os.path.abspath(os.path.join(os.path.dirname(__file__), 'regev_lab_cell_cycle_genes.txt')))]

    # Define two lists, genes associated to the S phase and genes associated to the G2M phase
    s_genes = cell_cycle_genes[:43]
    g2m_genes = cell_cycle_genes[43:]
    cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]

    # Perform cell cycle scoring. Cell cycle scoring adds three slots in data, a score for S phase, a score for G2M phase and the predicted cell cycle phase.
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

    adata_cc_genes = adata[:, cell_cycle_genes]
    sc.tl.pca(adata_cc_genes)
    # sc.pl.pca_scatter(adata_cc_genes, color='phase')

    # Regress out both S score and G2M score
    sc.pp.regress_out(adata, ['S_score', 'G2M_score'])
    sc.pp.scale(adata)

    # Reproject dataset using cell cycle genes again. Since we regressed the scores, no effect of cell cycle is now evident.
    adata_cc_genes = adata[:, cell_cycle_genes]
    sc.tl.pca(adata_cc_genes)
    # sc.pl.pca_scatter(adata_cc_genes, color='phase')

    return adata


def run_dimension_reduction(adata, layer=None, n_neighbors=10, n_pcs=40, resolution=1):
    if layer is not None:
        adata_temp = adata.copy()
        adata_temp.X = adata_temp.layers[layer]
        adata_temp = adata

        # Principal component analysis
        sc.tl.pca(adata_temp, svd_solver='arpack')
        adata.uns['pca'] = adata_temp.uns["pca"]
        adata.obsm[layer+'_pca'] = adata_temp.obsm["X_pca"]

        # Computing the neighborhood graph
        sc.pp.neighbors(adata_temp, n_neighbors=n_neighbors, n_pcs=n_pcs)
        adata.uns['neighbors'] = adata_temp.uns["neighbors"]
        adata.obsp['distances'] = adata_temp.obsp['distances']
        adata.obsp['connectivities'] = adata_temp.obsp['connectivities']

        # tSNE
        sc.tl.tsne(adata_temp)
        adata.uns['tsne'] = adata_temp.uns["tsne"]
        adata.obsm[layer+'_tsne'] = adata_temp.obsm["X_tsne"]

        # Clustering the neighborhood graph
        sc.tl.umap(adata_temp)
        adata.uns['umap'] = adata_temp.uns["umap"]
        adata.obsm[layer+'_umap'] = adata_temp.obsm["X_umap"]

        sc.tl.leiden(adata_temp, resolution=resolution)
        adata.uns['leiden'] = adata_temp.uns["leiden"]
        sc.tl.louvain(adata_temp, resolution=resolution)
        adata.uns['louvain'] = adata_temp.uns["louvain"]
        adata_temp = None
    else:
        # Principal component analysis
        sc.tl.pca(adata, svd_solver='arpack')

        # Computing the neighborhood graph
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

        # tSNE
        sc.tl.tsne(adata)

        # Clustering the neighborhood graph
        sc.tl.umap(adata)
        sc.tl.leiden(adata, resolution=resolution)
        sc.tl.louvain(adata, resolution=resolution)

    return adata