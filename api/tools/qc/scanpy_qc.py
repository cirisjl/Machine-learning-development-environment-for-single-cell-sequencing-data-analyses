import numpy as np
import pandas as pd
import scanpy as sc
import sklearn
import scipy
# sys.path.append('..')
from tools.formating.formating import *
sc.settings.verbosity=3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')


def scanpy_qc(adata):
        adata.var_names_make_unique() 

        # Preprocessing
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=3)
        adata.var['mt']=adata.var_names.str.startswith('MT-')  # Annotate the group of mitochondrial genes as 'mt'
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
        adata=adata[adata.obs.n_genes_by_counts < 2500, :]
        adata=adata[adata.obs.pct_counts_mt < 5, :]
        sc.pp.normalize_total(adata, target_sum=1e4)

        sc.pp.log1p(adata)
        
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)


        adata=adata[:, adata.var.highly_variable] # Do the filtering
        sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
        sc.pp.scale(adata, max_value=10)

        # Principal component analysis
        sc.tl.pca(adata, svd_solver='arpack')
        
        # Computing the neighborhood graph
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

        # Clustering the neighborhood graph
        sc.tl.umap(adata)
        sc.tl.leiden(adata)
        sc.tl.leiden(adata, resolution=1.5, key_added="cluster2")
        sc.pl.umap(adata, color=['leiden','cluster2'])

        # return adata, output
        return adata

