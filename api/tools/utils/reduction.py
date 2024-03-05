import warnings
import scanpy as sc
from anndata import AnnData



def run_dimension_reduction(adata, layer=None, n_neighbors=15, use_rep=None, n_pcs=None, resolution=1):
    if layer is not None:
        adata_temp = adata.copy()
        adata_temp.X = adata_temp.layers[layer]

        # Principal component analysis
        sc.tl.pca(adata_temp, svd_solver='arpack')
        adata.uns['pca'] = adata_temp.uns["pca"].copy()
        adata.obsm[layer+'_pca'] = adata_temp.obsm["X_pca"].copy()

        if use_rep is not None and adata_temp.obsm[use_rep].shape[1] < n_pcs:
            print(use_rep + " does not have enough Dimensions. Set n_pcs to " + adata_temp.obsm[use_rep].shape[1])
            n_pcs = adata_temp.obsm[use_rep].shape[1]

        # Computing the neighborhood graph
        sc.pp.neighbors(adata_temp, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep=use_rep)
        adata.uns['neighbors'] = adata_temp.uns["neighbors"].copy()
        adata.obsp['distances'] = adata_temp.obsp['distances'].copy()
        adata.obsp['connectivities'] = adata_temp.obsp['connectivities'].copy()

        # tSNE
        sc.tl.tsne(adata_temp)
        adata.uns['tsne'] = adata_temp.uns["tsne"].copy()
        adata.obsm[layer+'_tsne'] = adata_temp.obsm["X_tsne"].copy()

        # Clustering the neighborhood graph
        sc.tl.umap(adata_temp)
        adata.uns['umap'] = adata_temp.uns["umap"].copy()
        adata.obsm[layer+'_umap'] = adata_temp.obsm["X_umap"].copy()

        sc.tl.leiden(adata_temp, resolution=resolution)
        adata.uns[layer + '_leiden'] = adata_temp.uns["leiden"].copy()
        adata.obs[layer + '_leiden'] = adata_temp.obs["leiden"].copy()
        sc.tl.louvain(adata_temp, resolution=resolution)
        adata.uns[layer + '_louvain'] = adata_temp.uns["louvain"].copy()
        adata.obs[layer + '_louvain'] = adata_temp.obs["louvain"].copy()
        adata_temp = None
    else:
        # Principal component analysis
        sc.tl.pca(adata, svd_solver='arpack')

        # Computing the neighborhood graph
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep=use_rep)

        # tSNE
        sc.tl.tsne(adata)

        # Clustering the neighborhood graph
        sc.tl.umap(adata)
        leiden_key = "leiden"
        louvain_key = "louvain"
        if use_rep is not None:
            leiden_key = "leiden_" + use_rep
            louvain_key = "louvain_" + use_rep
        sc.tl.leiden(adata, key_added = leiden_key, resolution=resolution)
        sc.tl.louvain(adata, key_added = louvain_key, resolution=resolution)

    return adata
