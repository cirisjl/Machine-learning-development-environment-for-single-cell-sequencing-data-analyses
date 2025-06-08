import anndata
import scanpy as sc
import scprep
import sklearn.metrics


# test_data = adata.obsm["test"]
# denoised_data = adata.obsm["denoised"]
# train_data = adata.obsm["train"]
def imputation_metrics(adata, train, denoised, test):
    #Mean-squared error
    test_adata = anndata.AnnData(X=adata.layers[test], obs=adata.obs, var=adata.var)
    denoised_adata = anndata.AnnData(
        X=adata.layers[denoised], obs=adata.obs, var=adata.var
    )

    # scaling and transformation
    target_sum = 10000

    sc.pp.normalize_total(test_adata, target_sum)
    sc.pp.log1p(test_adata)

    sc.pp.normalize_total(denoised_adata, target_sum)
    sc.pp.log1p(denoised_adata)

    mse = sklearn.metrics.mean_squared_error(
        scprep.utils.toarray(test_adata.X), scprep.utils.toarray(denoised_adata.X)
    )

    # Poisson loss
    test_data = adata.layers[test]
    denoised_data = adata.layers[denoised]

    # scaling
    initial_sum = adata.layers[train].sum()
    target_sum = test_data.sum()
    denoised_data = denoised_data * target_sum / initial_sum

    possion = poisson_nll_loss(scprep.utils.toarray(test_data), scprep.utils.toarray(denoised_data))

    return mse, possion


def poisson_nll_loss(y_pred: np.ndarray, y_true: np.ndarray) -> float:
    return (y_pred - y_true * np.log(y_pred + 1e-6)).mean()