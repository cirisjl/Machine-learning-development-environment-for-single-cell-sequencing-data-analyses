import numpy as np
import pandas as pd
import scanpy as sc
import scrublet as scr

from tools.formating.formating import *

def predict_scrublet(path):
    adata = load_anndata(path)
    adata.var_names_make_unique()
    scrub = scr.Scrublet(adata.X, expected_doublet_rate = 0.076)
    adata.obs['doublet_scores'], adata.obs['predicted_doublets'] = scrub.scrub_doublets(min_counts=2, min_cells=3, 
                                                            min_gene_variability_pctl=85, n_prin_comps=30)
    
    scrub.plot_histogram()
    
    scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
    scrub.plot_embedding('UMAP', order_points=True)
    
    adata.obs['predicted_doublets'].value_counts()
    
    scrublet_path = os.path.join(os.path.dirname(path), "scrublet_calls.tsv")
    
    pd.DataFrame(adata.obs.iloc[:, -2:]).to_csv(scrublet_path,sep = '\t',header = False)
    
    return scrublet_path