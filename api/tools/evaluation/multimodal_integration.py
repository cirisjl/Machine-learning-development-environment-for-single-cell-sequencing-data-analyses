from scipy import sparse
import numpy as np
import scprep
import sklearn.decomposition
import sklearn.neighbors

# adata.obsm["aligned"]
# adata.obsm["mode2_aligned"]
def multimodal_integration_metrics(adata, aligned="aligned", mode2_aligned="mode2_aligned", proportion_neighbors=0.1, n_svd=100):
    # Mean squared error
    X = scprep.utils.toarray(adata.obsm["aligned"])
    Y = scprep.utils.toarray(adata.obsm["mode2_aligned"])

    X_shuffled = X[np.random.permutation(np.arange(X.shape[0])), :]
    error_random = np.mean(np.sum(_square(X_shuffled - Y)))
    error_abs = np.mean(np.sum(_square(X - Y)))
    mse = float('{:.4f}'.format(error_abs / error_random))

    # kNN Area Under the Curve
    n_svd = min([n_svd, min(adata.X.shape) - 1])
    n_neighbors = int(np.ceil(proportion_neighbors * adata.X.shape[0]))
    X_pca = sklearn.decomposition.TruncatedSVD(n_svd).fit_transform(adata.X)
    _, indices_true = (
        sklearn.neighbors.NearestNeighbors(n_neighbors=n_neighbors)
        .fit(X_pca)
        .kneighbors(X_pca)
    )
    _, indices_pred = (
        sklearn.neighbors.NearestNeighbors(n_neighbors=n_neighbors)
        .fit(adata.obsm[aligned])
        .kneighbors(adata.obsm[mode2_aligned])
    )
    neighbors_match = np.zeros(n_neighbors, dtype=int)
    for i in range(adata.shape[0]):
        _, pred_matches, true_matches = np.intersect1d(
            indices_pred[i], indices_true[i], return_indices=True
        )
        neighbors_match_idx = np.maximum(pred_matches, true_matches)
        neighbors_match += np.sum(
            np.arange(n_neighbors) >= neighbors_match_idx[:, None],
            axis=0,
        )

    neighbors_match_curve = neighbors_match / (
        np.arange(1, n_neighbors + 1) * adata.shape[0]
    )
    area_under_curve = float('{:.4f}'.format(np.mean(neighbors_match_curve)))

    return mse, area_under_curve


def _square(X):
    if sparse.issparse(X):
        X.data = X.data**2
        return X
    else:
        return scprep.utils.toarray(X) ** 2