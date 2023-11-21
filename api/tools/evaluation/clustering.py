from sklearn.metrics import adjusted_rand_score as ARI
from sklearn.metrics import normalized_mutual_info_score as NMI
from sklearn.metrics import silhouette_score

def clustering_scores(labels, labels_pred, embedding):
    asw_score = silhouette_score(embedding, labels)
    nmi_score = NMI(labels, labels_pred)
    ari_score = ARI(labels, labels_pred)
    print(
        "Clustering Scores:\nSilhouette: %.4f\nNMI: %.4f\nARI: %.4f"
        % (asw_score, nmi_score, ari_score)
    )
    return asw_score, nmi_score, ari_score