# scVI imports
import scvi
from tools.formating.formating import load_anndata, reset_x_to_raw
from exceptions.custom_exceptions import CeleryTaskException

def scvi_transfer(adata, refs = [], labels = None):
    if len(refs) > 0 and labels is not None:
        adata = reset_x_to_raw(adata)
        adata.obs['CellType'] = 'Unknown'
        adata.obs['Batch'] = 'Unknown'
        sc.pp.filter_genes(adata, min_cells = 10)
        adatas = [load_anndata(input) for input in refs]
        dater = sc.concat(adatas, join='outer')
        sc.pp.filter_genes(dater, min_cells = 10)
        dater = dater[~dater.obs[labels].isna()]
        dater = reset_x_to_raw(dater)
        dater.obs['CellType'] = dater.obs[labels]
        dater.obs['Batch'] = 'reference'
        dater = sc.concat((adata, dater))
        dater.obs['Sample'] = dater.obs.index

        sc.pp.highly_variable_genes(dater, flavor = 'seurat_v3', n_top_genes=3000, batch_key="Batch", subset = True)
        scvi.model.SCVI.setup_anndata(dater, batch_key='Batch', categorical_covariate_keys = ['Sample'])
        vae = scvi.model.SCVI(dater)
        vae.train(max_epochs = 400, early_stopping = True)

        lvae = scvi.model.SCANVI.from_scvi_model(vae, adata = dater, unlabeled_category = 'Unknown',
                                            labels_key = 'CellType')

        lvae.train(max_epochs=20, n_samples_per_label=100)

        dater.obs['scVI_predicted'] = lvae.predict(dater)
        dater.obs['scVI_transfer_score'] = lvae.predict(soft = True).max(axis = 1)
        dater = dater[dater.obs.Batch == 'Unknown']

        adata.obs = adata.obs.merge(right = dater.obs[['scVI_predicted', 'scVI_transfer_score']], left_index=True, right_index=True)
        adata.obs = adata.obs.drop('CellType', axis=1)
    else:
        raise CeleryTaskException(f"scVI annotation is failed due to empty user reference ({refs}) or cell labels ({labels}).")

    return adata