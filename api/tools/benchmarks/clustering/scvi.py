# Seed for reproducibility
import torch
import numpy as np
import pandas as pd
import api.tools.benchmarks.clustering.scanpy as sc
from typing import Tuple

import sys
sys.path.append('..')
# from tools.formating.formating import *
from tools.utils.utils import run_dimension_reduction
from tools.evaluation.monitor import *
from tools.evaluation.clustering import clustering_scores

# scVI imports
import scvi
from scvi.model.utils import mde
import pymde

torch.manual_seed(0)
np.random.seed(0)
sc.settings.verbosity = 0  # verbosity: errors (0), warnings (1), info (2), hints (3)


def clustering(adata, labels):
    # Start monitoring
    monitor = Monitor(1)

    model = scvi.model.SCVI(adata)
    model.train()

    latent = model.get_latent_representation()
    adata.obsm["X_scVI"] = latent
    denoised = model.get_normalized_expression(adata, library_size=1e4)
    adata.layers["scvi_normalized"] = model.get_normalized_expression(
        library_size=10e4
    )

    adata = run_dimension_reduction(adata, use_rep='X_scVI')

    adata.obsm["X_mde"] = mde(adata.obsm["X_scVI"])
    
    # Stop monitoring
    time_points, cpu_usage, mem_usage, gpu_mem_usage = monitor.stop()

    asw_score, nmi_score, ari_score = clustering_scores(adata.obs[labels], adata.obs["leiden_X_scVI"], adata.obsm["X_mde"])

    return asw_score, nmi_score, ari_score, time_points, cpu_usage, mem_usage, gpu_mem_usage



