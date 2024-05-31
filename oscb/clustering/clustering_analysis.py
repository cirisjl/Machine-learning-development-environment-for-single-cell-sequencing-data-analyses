# clustering.py

import scanpy as sc
import pandas as pd
import numpy as np
import argparse
from sklearn.metrics import adjusted_rand_score as ARI
from sklearn.metrics import normalized_mutual_info_score as NMI
from sklearn.metrics import silhouette_score

def clustering_scores(labels, labels_pred, embedding):
    if len(np.unique(labels)) < 2:
        print("Error: Number of unique labels is less than 2. Cannot compute silhouette score.")
        return None, None, None

    asw_score = silhouette_score(embedding, labels)
    nmi_score = NMI(labels, labels_pred)
    ari_score = ARI(labels, labels_pred)
    asw_score = float('{:.4f}'.format(asw_score))
    nmi_score = float('{:.4f}'.format(nmi_score))
    ari_score = float('{:.4f}'.format(ari_score))

    print(
        "Clustering Scores:\nSilhouette: %.4f\nNMI: %.4f\nARI: %.4f"
        % (asw_score, nmi_score, ari_score)
    )
    return asw_score, nmi_score, ari_score

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Process .h5ad file for clustering scores.')
    parser.add_argument('file_path', type=str, help='Path to the .h5ad file')
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Load the .h5ad file
    adata = sc.read_h5ad(args.file_path)

    # Extract the required columns, assuming 'cell_type', 'cell_states', and 'X_umap' are correct
    labels = adata.obs['cell_type']
    labels_pred = adata.obs['cell_states'] 
    embedding = adata.obsm['X_umap']

    # Print unique labels to check for issues
    print("Unique labels:", np.unique(labels))

    # Calculate clustering scores for the chosen clustering labels
    print("Clustering scores:")
    clustering_scores(labels.values, labels_pred.values, embedding)
