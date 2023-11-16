# import dash
# import dash_core_components as dcc
# import dash_html_components as html
# import plotly.graph_objs as go

# import matplotlib.pyplot as plt
# import seaborn as sns

import json 
# import scanpy as sc
import pandas as pd
import numpy as np
from anndata import AnnData

from tools.formating.plotConstants import *


def plot_UMAP(adata, clustering_plot_type="n_genes", selected_cell_intersection=[], n_dim=2): # clustering_plot_type: 'n_genes', 'cluster.ids', 'leiden', 'louvain', 'seurat_clusters'
    print("[DEBUG] generating new UMAP plot")

    obs = adata.obs
    obsm = adata.obsm

    # validate that there is a 3D projection available if that was requested
    if (("X_umap_3D" in obsm.keys()) and (n_dim == 3)):
        coords = pd.DataFrame(obsm["X_umap_3D"], index=obs.index)
    else:
        n_dim = 2
        coords = pd.DataFrame(obsm["X_umap"], index=obs.index)
    
    traces = []
    for i, val in enumerate(sorted(obs[clustering_plot_type].unique())):
        a = obs[obs[clustering_plot_type] == val]
        b = coords[obs[clustering_plot_type] == val]
        s = []
        if (selected_cell_intersection in [None, [], ""]):
            s = list(range(len(a.index)))
        else:
            for c in selected_cell_intersection:
                if (c in a.index):
                    s.append(a.index.get_loc(c))
        if (n_dim == 2):
            traces.append({
                "type": "scattergl",
                "x": b[0].tolist(),
                "y": b[1].tolist(),
                "text": ["Cell ID: " + str(cell_id) for cell_id in a.index.astype(str)],
                "selectedpoints": s,
                "mode":'markers',
                "marker": {
                    'size': point_size_2d,
                    'line': {'width': point_line_width_2d, 'color': 'grey'},
                    "color": discrete_colors_3[i % len(discrete_colors_3)]
                },
                "unselected": {
                    "marker": {"opacity": min_opacity}
                },
                "selected": {
                    "marker": {"opacity": max_opacity}
                },
                "name": f"Cluster {val}"
            })
        elif (n_dim == 3):
            traces.append({
                "type": "scattergl",
                "x": b[0].tolist(),
                "y": b[1].tolist(),
                "z": b[2].tolist(),
                "text": ["Cell ID: " + str(cell_id) for cell_id in a.index.astype(str)],
                "selectedpoints": s,
                "mode":'markers',
                "marker": {
                    'size': point_size_3d,
                    'line': {'width': point_line_width_3d, 'color': 'grey'},
                    "color": discrete_colors_3[i % len(discrete_colors_3)]
                },
                "name": f"Cluster {val}"
            })
    if (n_dim == 2):
        return json.dumps({
            'data': traces,
            'layout': dict(
                xaxis={"title": "UMAP 1"},
                yaxis={"title": "UMAP 2"},
                margin=margin,
                # legend={'x': 0, 'y': 1},
                hovermode='closest',
                transition = {'duration': 250},
                autosize=True
                #width=4 * scale,
                #height=3 * scale
            )
        })
    elif (n_dim == 3):
        return json.dumps({
            'data': traces,
            'layout': dict(
                xaxis={"title": "UMAP 1"},
                yaxis={"title": "UMAP 2"},
                zaxis={"title": "UMAP 3"},
                margin=margin,
                # legend={'x': 0, 'y': 1},
                hovermode='closest',
                transition = {'duration': 250},
                autosize=True
                #width=4 * scale,
                #height=3 * scale
            )
        })


def plot_violin(adata, features=['n_counts', 'n_genes', 'pct_counts_mt', 'pct_counts_rb'], show_points=False):
    var    = adata.var
    obs    = adata.obs

    traces = []
    
    x_pos  = 1

    # n_traces = len(features)
    
    for i in features:
        if not i in obs.columns:
            print("[DEBUG] feature " + str(i) + " not in obs columns; skipping")
            continue
        if (show_points == False):
            traces.append({
                "type": "violin",
                "y": adata.obs_vector(i).tolist(), 
                "text": ["Cell ID: " + str(cell_id) for cell_id in obs.index.astype(str)],
                "opacity": 0.7,
                "box": {
                    "visible": True
                },
                "meanline": {
                    "visible": True
                },
                "points": "none",
                "name": str(i)
            })
            x_pos += 1
        
        elif (show_points == "all"):
            #kernel = gaussian_kde(adata.obs_vector(i))
            jittered_x = x_pos + 0.1 * np.random.standard_normal(len(adata.obs_vector(i)))

            traces.append({
                "type": "violin",
                "x": jittered_x.tolist(),
                "y": adata.obs_vector(i).tolist(), 
                "text": ["Cell ID: " + str(cell_id) for cell_id in obs.index.astype(str)],
                "mode": "markers",
                "opacity": 0.7,
                "marker": {
                    'size': point_size_2d
                },
                "name": str(i)
            })
            x_pos += 1

    if (traces in [[], None, ""]):
        print("[DEBUG] no traces added to violin plot")

    return json.dumps({
        'data': traces,
        'layout': dict(
            # xaxis={"title": "Selected metadata features"},
            # yaxis={"title": "value"},
            margin=margin,
            # legend={'x': 0, 'y': 1},
            hovermode='closest',
            transition = {'duration': 100},
            autosize=True
            #width=4 * scale,
            #height=3 * scale
        )
    })


def plot_scatter(adata, feature1 = "n_counts", feature2 = "n_genes"):
    obs = adata.obs

    traces = []

    if feature1 in obs.columns and feature2 in obs.columns:
        traces.append({
                "type": "scatter",
                "x": adata.obs_vector(feature1).tolist(), 
                "y": adata.obs_vector(feature2).tolist(), 
                "text": ["Cell ID: " + str(cell_id) for cell_id in obs.index.astype(str)],
                "mode":'markers',
                "marker": {
                    'size': point_size_2d,
                    'line': {'width': point_line_width_2d, 'color': 'grey'}
                }
            })  
    else:
        print("[DEBUG] " + feature1 + " or " + feature2 + " is not in adata.obs.columns")

    return json.dumps({
        'data': traces,
        'layout': dict(
            xaxis={"title": feature1},
            yaxis={"title": feature2},
            margin=margin,
            # legend={'x': 0, 'y': 1},
            hovermode='closest',
            transition = {'duration': 100},
            autosize=True
            #width=4 * scale,
            #height=3 * scale
        )
    })


def plot_highest_expr_genes(adata, n_top=30):
    import scanpy as sc
    from scipy.sparse import issparse

    # compute the percentage of each gene per cell
    norm_dict = sc.pp.normalize_total(adata, target_sum=100, inplace=False)

    # identify the genes with the highest mean
    if issparse(norm_dict["X"]):
        mean_percent = norm_dict["X"].mean(axis=0).A1
        top_idx = np.argsort(mean_percent)[::-1][:n_top]
        counts_top_genes = norm_dict["X"][:, top_idx].A
    else:
        mean_percent = norm_dict["X"].mean(axis=0)
        top_idx = np.argsort(mean_percent)[::-1][:n_top]
        counts_top_genes = norm_dict["X"][:, top_idx]
    columns = (
        adata.var_names[top_idx]
    )
    counts_top_genes = pd.DataFrame(
        counts_top_genes, index=adata.obs_names, columns=columns
    )

    y = counts_top_genes.columns.tolist()

    traces = []

    for i in y:
        traces.append({
                "type": "box",
                "x": counts_top_genes[i].tolist(), 
                "y": i, 
                "boxpoints": "all",
                # "jitter": 0.5,
                # "whiskerwidth": 0.2,
                "marker": {
                    "size": 2,
                    "line": {"width": 1, "color": 'grey'}
                },
                "name": str(i)
            })  

    return json.dumps({
        'data': traces,
        'layout': dict(
            title='% of total counts',
            xaxis={
                "autorange": True,
                "showgrid": True,
                "zeroline": True,
                "gridcolor": 'grey',
                "gridwidth": 1,
                "zerolinecolor": 'grey',
                "zerolinewidth": 2
            },
            margin=margin,
            # legend={'x': 0, 'y': 1},
            # hovermode='closest',
            transition = {'duration': 100},
            autosize=True,
            showlegend=False,
            width=4*scale,
            height=3*scale
        )
    })