import os
import muon
import numpy as np
import scanpy as sc
from scipy.sparse import csr_matrix
import scvi
import seaborn as sns
import torch
import scipy
from scipy import sparse
import sys
sys.path.append('..')
from tools.formating.formating import get_scvi_path

scvi.settings.seed = 0
torch.set_float32_matmul_precision("high")

def run_multivi(mdata_path, rna_subset="rna_subset", atac_subset="atac_subset"):
    mdata = muon.read_h5mu(mdata_path)
    model_dir = get_scvi_path(mdata_path, "multivi")
    model = None
    if not os.path.exists(model_dir):
        scvi.model.MULTIVI.setup_mudata(
            mdata,
            modalities={
                "rna_layer": rna_subset,
                "atac_layer": atac_subset,
            },
        )

        model = scvi.model.MULTIVI(
            mdata,
            n_genes=len(mdata.mod[rna_subset].var),
            n_regions=len(mdata.mod[atac_subset].var),
        )

        # For our sparse matrices, we want CSR rather than CSC as training will be faster
        if type(mdata.mod[rna_subset].X) == scipy.sparse._csc.csc_matrix:
            mdata.mod[rna_subset].X = mdata.mod[rna_subset].X.tocsr()
        elif type(mdata.mod[rna_subset].X) == np.matrix:
            mdata.mod[rna_subset].X = sparse.csr_matrix(mdata.mod[rna_subset].X)
        if type(mdata.mod[atac_subset].X) == scipy.sparse._csc.csc_matrix:
            mdata.mod[atac_subset].X = mdata.mod[atac_subset].X.tocsr()
        elif type(mdata.mod[atac_subset].X) == np.matrix:
            mdata.mod[atac_subset].X = sparse.csr_matrix(mdata.mod[atac_subset].X)
        mdata.update()
        model.train()
        model.save(model_dir, overwrite=True)
    else: 
        model = scvi.model.MULTIVI.load(model_dir, adata=mdata)

    mdata = model.adata

    # Below we an cell annotations for modality, so we can color the UMAP

    MULTIVI_LATENT_KEY = "X_multivi"

    mdata.obsm[MULTIVI_LATENT_KEY] = model.get_latent_representation()
    sc.pp.neighbors(mdata, use_rep=MULTIVI_LATENT_KEY)
    sc.tl.umap(mdata, min_dist=0.2)

    n = mdata.n_obs // 3

    # initialize the column first
    mdata.obs["modality"] = ""

    # set modality of first third to rna
    mdata.obs.iloc[:n, mdata.obs.columns.get_loc("modality")] = "expression"

    # set modality of second third to both
    mdata.obs.iloc[n : 2 * n, mdata.obs.columns.get_loc("modality")] = "paired"

    # set modality of last third to atac
    mdata.obs.iloc[2 * n :, mdata.obs.columns.get_loc("modality")] = "accessibility"

    sc.pl.umap(mdata, color="modality")
    # Impute missing modalities
    mdata.obsm["MultiVI_imputed"] = model.get_normalized_expression()

    return mdata