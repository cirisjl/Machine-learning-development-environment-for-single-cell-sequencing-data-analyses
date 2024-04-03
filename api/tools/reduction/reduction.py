import warnings
import scanpy as sc
from anndata import AnnData


def run_dimension_reduction(adata, layer=None, n_neighbors=15, use_rep=None, n_pcs=None, random_state=0):
    msg = None
    if layer == "Pearson_residuals":
        msg = "Normalize Pearson_residuals may create NaN values, which are not accepted by PCA."
        return adata, msg

    if layer is not None and layer+'_umap' not in adata.obsm.keys():
        adata_temp = adata.copy()
        adata_temp.X = adata_temp.layers[layer]

        # Principal component analysis
        sc.tl.pca(adata_temp, svd_solver='arpack', random_state=random_state)
        adata.uns['pca'] = adata_temp.uns["pca"].copy()
        adata.obsm[layer+'_pca'] = adata_temp.obsm["X_pca"].copy()

        if use_rep is not None and adata_temp.obsm[use_rep].shape[1] < n_pcs:
            msg = f"{use_rep} does not have enough Dimensions. Set n_pcs to {adata_temp.obsm[use_rep].shape[1]}."
            n_pcs = adata_temp.obsm[use_rep].shape[1]

        # Computing the neighborhood graph
        sc.pp.neighbors(adata_temp, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep=use_rep, random_state=random_state)
        adata.uns['neighbors'] = adata_temp.uns["neighbors"].copy()
        adata.obsp['distances'] = adata_temp.obsp['distances'].copy()
        adata.obsp['connectivities'] = adata_temp.obsp['connectivities'].copy()

        # tSNE
        sc.tl.tsne(adata_temp)
        adata.uns['tsne'] = adata_temp.uns["tsne"].copy()
        adata.obsm[layer+'_tsne'] = adata_temp.obsm["X_tsne"].copy()

        # 2D UMAP
        sc.tl.umap(adata_temp, random_state=random_state, 
                    init_pos="spectral", n_components=2, 
                    copy=False, maxiter=None)
        adata.uns['umap'] = adata_temp.uns["umap"].copy()
        adata.obsm[layer+'_umap'] = adata_temp.obsm["X_umap"].copy()

        # 3D UMAP
        adata_3D = sc.tl.umap(adata_temp, random_state=random_state, 
                          init_pos="spectral", n_components=3, 
                          copy=True, maxiter=None)
        adata.obsm[layer+"_umap_3D"] = adata_3D.obsm["X_umap"]
        adata_3D = None

    elif 'X_umap' not in adata.obsm.keys():
        # Principal component analysis
        sc.tl.pca(adata, svd_solver='arpack', random_state=random_state)

        # Computing the neighborhood graph
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep=use_rep, random_state=random_state)

        # tSNE
        sc.tl.tsne(adata)

        # 2D UMAP
        sc.tl.umap(adata, random_state=random_state, 
                    init_pos="spectral", n_components=2, 
                    copy=False, maxiter=None)
        
        # 3D UMAP
        adata_3D = sc.tl.umap(adata, random_state=random_state, 
                          init_pos="spectral", n_components=3, 
                          copy=True, maxiter=None)
        adata.obsm["X_umap_3D"] = adata_3D.obsm["X_umap"]
        adata_3D = None
    else:
        msg = "{layer}_umap already exists, skipped."

    return adata, msg


def run_clustering(adata, layer=None, use_rep=None, resolution=1, random_state=0):
    if layer == "Pearson_residuals":
        print("Normalize Pearson_residuals may create NaN values, which are not accepted by PCA.")
        return adata
    
    if layer is not None and layer + '_louvain' not in adata.obs.keys():
        # Clustering the neighborhood graph
        adata_temp = adata.copy()
        sc.tl.leiden(adata_temp, resolution=resolution, 
                    random_state=random_state, n_iterations=3)
        adata.uns[layer + '_leiden'] = adata_temp.uns["leiden"].copy()
        adata.obs[layer + '_leiden'] = adata_temp.obs["leiden"].copy()
        sc.tl.louvain(adata_temp, resolution=resolution)
        adata.uns[layer + '_louvain'] = adata_temp.uns["louvain"].copy()
        adata.obs[layer + '_louvain'] = adata_temp.obs["louvain"].copy()
        adata_temp = None
    elif 'louvain' not in adata.obs.keys():
        # Clustering the neighborhood graph
        leiden_key = "leiden"
        louvain_key = "louvain"
        if use_rep is not None:
            leiden_key = "leiden_" + use_rep
            louvain_key = "louvain_" + use_rep
        sc.tl.leiden(adata, key_added = leiden_key, resolution=resolution, 
                    random_state=random_state, n_iterations=3)
        sc.tl.louvain(adata, key_added = louvain_key, resolution=resolution)
    else:
        print("Cluster for {layer} already exists, skipped.")

    return adata
