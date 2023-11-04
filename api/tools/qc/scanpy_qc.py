import numpy as np
import pandas as pd
import scanpy as sc
import warnings
warnings.filterwarnings('ignore')
import sklearn
import scipy
# sys.path.append('..')
from tools.formating.formating import *
sc.settings.verbosity=3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')


def run_scanpy_qc(adata, min_genes=200, min_cells=3, target_sum=1e4, regress_cell_cycle=False):
        if adata is None:
            print("File format is not supported.")
            return None

        adata.var_names_make_unique() 

        # Basic filtering
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=min_cells)
        adata.var['mt']=adata.var_names.str.startswith('MT-')  # Annotate the group of mitochondrial genes as 'mt'
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
        adata=adata[adata.obs.n_genes_by_counts < 2500, :]
        adata=adata[adata.obs.pct_counts_mt < 5, :]
        sc.pp.normalize_total(adata, target_sum=target_sum)

        sc.pp.log1p(adata)

        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)


        adata=adata[:, adata.var.highly_variable] # Do the filtering
        sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
        sc.pp.scale(adata, max_value=10)

        # Regress both S score and G2M score for cell cycle
        if(regress_cell_cycle):
             adata = regress_cell_cycle(adata)

        # Principal component analysis
        sc.tl.pca(adata, svd_solver='arpack')

        # Computing the neighborhood graph
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

        # tSNE
        sc.tl.tsne(adata)

        # Clustering the neighborhood graph
        sc.tl.umap(adata)
        sc.tl.leiden(adata)
        sc.tl.leiden(adata, resolution=1.5, key_added="cluster2")
        sc.pl.umap(adata, color=['leiden','cluster2'])

        # return adata, output
        return adata


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