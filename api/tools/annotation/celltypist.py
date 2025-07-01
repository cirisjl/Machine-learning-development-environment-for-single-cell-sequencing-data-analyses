import celltypist
import scanpy as sc
import os
from celltypist import models
from tools.formating.formating import load_anndata, reset_x_to_raw
models.get_all_models()


def run_celltypist(adata, model_name, refs = [], labels = None, species = 'mouse'):
    model = celltypist.Model.load(model_name)
    if species == 'mouse' and "Mouse" not in model_name:
        model.convert()
    adata = reset_x_to_raw(adata)

    sc.pp.filter_genes(adata, min_cells = 10)
    sc.pp.normalize_total(adata, target_sum=1e4) #not recommended for typical pp
    sc.pp.log1p(adata)
    
    if type(adata.X) != np.ndarray:
        adata.X = adata.X.toarray()
    
    predictions = celltypist.annotate(adata, model=model, majority_voting=True)
    predictions_adata = predictions.to_adata()
    adata.obs["celltypist_label"] = predictions_adata.obs.loc[adata.obs.index, "predicted_labels"]
    adata.obs["celltypist_score"] = predictions_adata.obs.loc[adata.obs.index, "conf_score"]

    if len(refs) > 0 and labels is not None:
        for input in refs:
            try:
                name = os.path.basename(input).split(".")[0]
                ref_ad = load_anndata(input)
                ref_ad = reset_x_to_raw(ref_ad)
                sc.pp.filter_genes(ref_ad, min_cells = 10)
                sc.pp.normalize_total(ref_ad, target_sum = 1e4) #Note this is only for cell annotation, recommended by authors but not best
                sc.pp.log1p(ref_ad)

                ref_ad = ref_ad[~ref_ad.obs[labels].isna()]
                ref_model = celltypist.train(ref_ad, labels = labels, n_jobs = 4, use_SGD = False, feature_selection = True, top_genes = 300)
                ref_predictions = celltypist.annotate(adata, model=ref_model, majority_voting=False)
                ref_predictions_adata = ref_predictions.to_adata()
                adata.obs["ref_"+name+"_label"] = ref_predictions_adata.obs.loc[adata.obs.index, "predicted_labels"]
                adata.obs["ref_"+name+"_score"] = ref_predictions_adata.obs.loc[adata.obs.index, "conf_score"]
            except Exception as e:
                print(e)
                continue

    return adata
            
    
