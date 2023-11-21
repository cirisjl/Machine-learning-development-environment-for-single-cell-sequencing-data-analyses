""" 
Calculate evaluation score from to get the text of file

return ARI, homogeneity, and sil
example run: calc_cluster_benchmark(adata, label_name_predict = 'leiden', label_name_benchmark = 'cluster2')
"""
def calc_cluster_benchmark(adata, label_name_predict, label_name_benchmark, reduction = 'X_pca'):
    ari_score = sklearn.metrics.adjusted_rand_score(adata.obs[label_name_predict], adata.obs[label_name_benchmark])
    homogeneity_score = sklearn.metrics.homogeneity_completeness_v_measure(adata.obs[label_name_predict], adata.obs[label_name_benchmark])[2]
    sil=sklearn.metrics.silhouette_score(adata.obsm[reduction], adata.obs[label_name_predict])
    return [ari_score, homogeneity_score, sil]