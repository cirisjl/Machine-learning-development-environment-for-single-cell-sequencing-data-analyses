import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go

import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
import pandas as pd
import numpy as np
from anndata import AnnData

from tools.formating.plotConstants import *


def plot_UMAP(adata, clustering_plot_type="cell_type", selected_cell_intersection=[], n_dim=2):
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
        if (selected_cell_intersection in [None, []]):
            s = list(range(len(a.index)))
        else:
            for c in selected_cell_intersection:
                if (c in a.index):
                    s.append(a.index.get_loc(c))
        if (n_dim == 2):
            traces.append(
                go.Scattergl(
                    x=b[0],
                    y=b[1],
                    text="Cell ID: " + a.index,
                    mode='markers',
                    selectedpoints=s,
                    marker={
                        'size': point_size_2d,
                        'line': {'width': point_line_width_2d, 'color': 'grey'},
                        "color": discrete_colors_3[i%len(discrete_colors_3)]
                    },
                    unselected={
                        "marker": {"opacity": min_opacity,
                        }
                    },
                    selected={
                        "marker": {"opacity": max_opacity,
                        }
                    },
                    name=("Cluster " + str(val))
                )
            )
        elif (n_dim == 3):
            traces.append(
                go.Scatter3d(
                    x=b[0],
                    y=b[1],
                    z=b[2],
                    text="Cell ID: " + a.index,
                    mode='markers',
                    marker={
                        'size': point_size_3d,
                        'line': {'width': point_line_width_3d, 'color': 'grey'},
                        "color": discrete_colors_3[i%len(discrete_colors_3)]
                    },
                    name=("Cluster " + str(val))
                )
            )
    if (n_dim == 2):
        return {
            'data': traces,
            'layout': dict(
                xaxis={"title": "UMAP 1"},
                yaxis={"title": "UMAP 2"},
                margin=margin,
                legend={'x': 0, 'y': 1},
                hovermode='closest',
                transition = {'duration': 250},
                autosize=True
                #width=4 * scale,
                #height=3 * scale
            )
        }
    elif (n_dim == 3):
        return {
            'data': traces,
            'layout': dict(
                xaxis={"title": "UMAP 1"},
                yaxis={"title": "UMAP 2"},
                zaxis={"title": "UMAP 3"},
                margin=margin,
                legend={'x': 0, 'y': 1},
                hovermode='closest',
                transition = {'duration': 250},
                autosize=True
                #width=4 * scale,
                #height=3 * scale
            )
        }


def plot_violin(adata, selected_genes=[], show_points = "all"):
    var    = adata.var
    obs    = adata.obs

    traces = []
    
    x_pos  = 1
    
    if (selected_genes in [None, []]):
        selected_genes = adata.var.index

    n_traces = len(selected_genes)
    
    for i in selected_genes:
        if not ((i in var.index) or (i in obs)):
            print("[DEBUG] gene " + str(i) + " not in var index; skipping")
            continue
        if (show_points == False):
            traces.append(
                go.Violin(
                    y=adata.obs_vector(i),
                    text="Cell ID: " + obs.index,
                    opacity=0.7,
                    box_visible=True,
                    meanline_visible=True,
                    points=False,
                    name=(str(i))
                )
            )
            x_pos += 1
        
        elif (show_points == "all"):
            #kernel = gaussian_kde(adata.obs_vector(i))
            jittered_x = x_pos + 0.1 * np.random.standard_normal(len(adata.obs_vector(i)))

            traces.append(
                go.Scattergl(
                    x=jittered_x,
                    y=adata.obs_vector(i),
                    text="Cell ID: " + obs.index,
                    mode="markers",
                    opacity=0.7,
                    marker={
                        'size': point_size_2d,
                    },
                    name=(str(i)),
                )
            )
            x_pos += 1

    if (traces in [[], None]):
        print("[DEBUG] no traces added to violin plot")
        return dash.no_update

    return {
        'data': traces,
        'layout': dict(
            xaxis={"title": "Gene/factor"},
            yaxis={"title": "Expression"},
            margin=margin,
            legend={'x': 0, 'y': 1},
            hovermode='closest',
            transition = {'duration': 100},
            autosize=True
            #width=4 * scale,
            #height=3 * scale
        )
    }