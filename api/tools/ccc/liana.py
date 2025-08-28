import numpy as np
import liana as li
import scanpy as sc
from tools.formating.formating import is_normalized
from liana.mt import rank_aggregate
from liana.method import singlecellsignalr, connectome, cellphonedb, natmi, logfc, cellchat, geometric_mean

# https://liana-py.readthedocs.io/en/latest/notebooks/basic_usage.html
def run_liana_ccc(adata, cell_type_label, species):
    if adata is None:
        raise ValueError("Failed to load AnnData object.")
    
    if is_normalized(adata.X, 200) and not check_nonnegative_integers(adata.X):
        redislogger.info(unique_id, "adata.X is not raw counts.")
        if adata.raw.X is not None:
            redislogger.info(unique_id, "Use adata.raw.X instead. Copy adata.X to layer 'normalized_X'.")
            adata.layers["normalized_X"] = adata.X.copy()
            adata.X = adata.raw.X.copy()
        elif "raw_counts" in adata.layers.keys():
            redislogger.info(unique_id, "Use layer 'raw_counts' instead. Copy adata.X to layer 'normalized_X'.")
            adata.layers["normalized_X"] = adata.X.copy()
            adata.X = adata.layers['raw_counts'].copy()
        else:
            raise ValueError("Liana only take raw counts, not normalized data.")

    resource_name = 'consensus'
    if species == "mouse":
        resource_name = "mouseconsensus"
            
    # run cellphonedb
    # cellphonedb(adata,
    #         groupby=cell_type_label, 
    #         # NOTE by default the resource uses HUMAN gene symbols
    #         resource_name='consensus',
    #         expr_prop=0.1,
    #         verbose=True, key_added='cpdb_res')

    # by default, liana's output is saved in place:
    # adata.uns['cpdb_res'].head()

    # Run rank_aggregate
    li.mt.rank_aggregate(adata, 
                        groupby=cell_type_label,
                        resource_name=resource_name,
                        expr_prop=0.1,
                        verbose=True)
    # adata.uns['liana_res'].head()

    # methods = [logfc, geometric_mean]
    # new_rank_aggregate = li.mt.AggregateClass(li.mt.aggregate_meta, methods=methods)
    # new_rank_aggregate(adata,
    #                 groupby='bulk_labels',
    #                 expr_prop=0.1, 
    #                 verbose=True,
    #                 # Note that with this option, we don't perform permutations
    #                 # and hence we exclude the p-value for geometric_mean, as well as specificity_rank
    #                 n_perms=None,
    #                 use_raw=True,
    #                 )
    
    # adata.uns['liana_res'].head()

    return adata