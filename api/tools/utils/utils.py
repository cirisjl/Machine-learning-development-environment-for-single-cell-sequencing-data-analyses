import warnings
import scanpy as sc
from typing import Union, Optional, Tuple, Collection, Sequence, Iterable, Literal
import numpy as np
from numpy import random
import scipy as sp
from scipy.sparse import issparse, isspmatrix_csr, csr_matrix, spmatrix
from anndata import AnnData

AnyRandom = Union[None, int, random.RandomState]  # random.Generator


def sc_train_test_split(
    data: Union[AnnData, np.ndarray, spmatrix],
    train_fraction: [float] = 0.8,
    n_obs: Optional[int] = None,
    random_state: AnyRandom = 0
    ) -> Optional[AnnData]:

    random.seed(random_state)
    total_obs = data.n_obs if isinstance(data, AnnData) else data.shape[0]
    obs_indices = np.arange(total_obs)

    if n_obs is not None:
        train_n_obs = n_obs
    elif train_fraction is not None:
        if train_fraction > 1 or train_fraction < 0:
            raise ValueError(f"`train_fraction` needs to be within [0, 1], not {train_fraction}")
        train_n_obs = int(train_fraction * total_obs)
    else:
        raise ValueError("Either pass `n_obs` or `train_fraction`.")

    train_indices = random.choice(total_obs, size=train_n_obs, replace=False)
    test_indices = np.setdiff1d(obs_indices, train_indices)

    if isinstance(data, AnnData):
        if data.isbacked:
            return data[train_indices].to_memory(), data[test_indices].to_memory()
        else:
            return data[train_indices].copy(), data[test_indices].copy()
    else:
        X = data
        return X[train_indices], X[test_indices]


def sc_train_val_test_split(
    data: Union[AnnData, np.ndarray, spmatrix],
    train_fraction: [float] = 0.8,
    validation_fraction: [float] = 0.1,
    test_fraction: [float] = 0.1
    ) -> Optional[AnnData]:

    if train_fraction is None or validation_fraction is None or test_fraction is None:
        raise ValueError("`train_fraction`, `validation_fraction` and `test_fraction` are required.")
    if train_fraction > 1 or train_fraction < 0:
        raise ValueError(f"`train_fraction` needs to be within [0, 1], not {train_fraction}")
    if validation_fraction > 1 or validation_fraction < 0:
        raise ValueError(f"`validation_fraction` needs to be within [0, 1], not {validation_fraction}")
    if test_fraction > 1 or test_fraction < 0:
        raise ValueError(f"`test_fraction` needs to be within [0, 1], not {test_fraction}")
    if train_fraction + validation_fraction + test_fraction != 1:
        raise ValueError("`train_fraction`, `validation_fraction` and `test_fraction` must be added up to 1.")

    train, validation_test = sc_train_test_split(data, train_fraction=train_fraction)
    validation, test = sc_train_test_split(validation_test, train_fraction=validation_fraction/(validation_fraction+test_fraction))

    return train, validation, test


def run_dimension_reduction(adata, layer=None, n_neighbors=10, use_rep=None, n_pcs=40, resolution=1):
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
        sc.tl.louvain(adata_temp, resolution=resolution)
        adata.uns[layer + '_louvain'] = adata_temp.uns["louvain"].copy()
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
        sc.tl.leiden(adata, resolution=resolution)
        sc.tl.louvain(adata, resolution=resolution)

    return adata
