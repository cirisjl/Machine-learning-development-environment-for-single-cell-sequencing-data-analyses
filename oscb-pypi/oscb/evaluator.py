from .evaluation.annotation import *
from .evaluation.ccc import *
from .evaluation.clustering import *
from .evaluation.imputation import *
from .evaluation.integration import *
from .evaluation.multimodal import *
from .evaluation.trajectory import *
from .evaluation.annotation import *
from datetime import datetime
from .utils import *
import requests
import json


def eval(adata, adata_int=None, benchmarks_id=None, task=None, cluster_key=None, label_key=None, label_pred_key=None, embedding_key=None, ccc_pred="ccc_pred", ccc_target="ccc_target", score="score", denoised_layer=None, train='train', test='test', mod1_key='rna', mod2_key='atac', traj_key=None, bm_traj_key=None, root_node=None, species=None, server_endpoint=server_endpoint+'benchmarks/', method="Your method"):
    if adata is None:
        raise ValueError("adata is required.")
    
    benchmarks = None
    current_date_and_time = datetime.now()
    benchmarks_data = None

    if benchmarks_id is not None:
        dataset_id, task = get_dataset_id(benchmarks_id)
        url = server_endpoint + benchmarks_id
        response = requests.get(url)
        if response.status_code == 200:
            try:
                benchmarks = response.json()
                benchmarks_data = benchmarks['benchmarks_plot']['data']
                match task:
                    case "Clustering" | "CL":
                        label_key = benchmarks['label']

                    case "Imputation" | "IM":
                        species = benchmarks['species']

                    case "Batch Integration" | "BI":
                        label_key = benchmarks['label']
                        batch_key = benchmarks['batch_key']
                        species = benchmarks['species']

                    case "Trajectory" | "TJ":
                        label_key = benchmarks['label']
                        root_node = benchmarks['origin_group']
                        bm_traj_key = benchmarks['bm_traj']

                    case "Cell-Cell Communication" | "CCC":
                        label_key = benchmarks['label']
                        ccc_target = benchmarks['ccc_target']
                        species = benchmarks['species']

                    case "Multimodal Data Integration" | "MI":
                        mod1_key = benchmarks['mod1']
                        mod2_key = benchmarks['mod2']
                        label_key = benchmarks['label']
                        batch_key = benchmarks['batch_key']

                    case "Cell Type Annotation" | "CT":
                        label_key = benchmarks['label']
                        # species = benchmarks['species']

            except Exception as e:
                print(f"Failed to get Benchmarks: {str(e)}")
        else:
            print(f"Failed to get Benchmarks: {benchmarks_id}.")

    if task is not None:
        task_info = {
            "benchmarksId": benchmarks_id,
            "datasetId": dataset_id,
            "task_type": task,
            "tool": method,
            "created_on": current_date_and_time
        }
        match task:
            case "Clustering" | "CL":
                if cluster_key is not None and label_key is not None and embedding_key is not None:
                    asw_score, nmi_score, ari_score, fm_score = clustering_metrics(adata.obs[label_key], adata.obs[cluster_key], adata.obsm[embedding_key])
                    results = {
                        "benchmarksId": benchmarks_id,
                        "datasetId": dataset_id,
                        "task_type": task,
                        "tool": method,
                        "Silhouette": asw_score,
                        "NMI": nmi_score,
                        "ARI": ari_score,
                        "Fowlkes Mallows": fm_score,
                        "created_on": current_date_and_time
                    }
                    if benchmarks_data is not None:
                        labels, y_labels, data = get_bar_plot_data(benchmarks_data, user_results=results)
                        plot_bars(task, labels, y_labels, data)

                    return results
                else: 
                    raise ValueError(f"cluster_key, label_key and embedding_key are required for {task}.")

            case "Imputation" | "IM":
                if denoised_layer is not None:
                    mse, possion = imputation_metrics(adata, denoised_layer=denoised_layer)
                    results = {
                        "benchmarksId": benchmarks_id,
                        "datasetId": dataset_id,
                        "task_type": task,
                        "tool": method,
                        "MSE": mse,
                        "Possion": possion,
                        "created_on": current_date_and_time
                    }
                    if benchmarks_data is not None:
                        labels, y_labels, data = get_bar_plot_data(benchmarks_data, user_results=results)
                        plot_bars(task, labels, y_labels, data)
                    return results
                else: 
                    raise ValueError(f"denoised_layer is required for {task}.")

            case "Batch Integration" | "BI":
                if adata_int is not None and label_key is not None and batch_key is not None:
                    metrics_dict = integration_metrics(adata, adata_int, batch_key=batch_key, label_key=label_key, species=species)
                    results = {**task_info, **metrics_dict}
                    if benchmarks_data is not None:
                        labels, y_labels, data = get_bar_plot_data(benchmarks_data, user_results=results)
                        plot_bars(task, labels, y_labels, data)
                    return results
                else: 
                    raise ValueError(f"adata_int, label_key and batch_key are required for {task}.")
                    
            case "Trajectory" | "TJ":
                if traj_key is not None and bm_traj_key is not None and root_node is not None:
                    ged_score, gks_score, jsc_score, ted_score, mean = trajectory_metrics(adata.uns[traj_key], adata.uns[bm_traj_key], adata.uns[root_node])
                    results = {
                        "benchmarksId": benchmarks_id,
                        "datasetId": dataset_id,
                        "task_type": task,
                        "tool": method,
                        "Graph Edit Distance": ged_score,
                        "Graph Kernel Score": gks_score,
                        "Jaccard Similarity Coefficient": jsc_score,
                        "Tree Edit Distance": ted_score,
                        "Mean": mean,
                        "created_on": current_date_and_time
                    }
                    if benchmarks_data is not None:
                        labels, y_labels, data = get_bar_plot_data(benchmarks_data, user_results=results)
                        plot_bars(task, labels, y_labels, data)
                    return results
                else: 
                    raise ValueError(f"adata_int, label_key and batch_key are required for {task}.")

            case "Cell-Cell Communication" | "CCC":
                if ccc_pred is not None and ccc_target is not None and score is not None:
                    auc_score, oddsratio_score = ccc_metrics(adata, ccc_pred=ccc_pred, ccc_target=ccc_target, score='score')
                    results = {
                        "benchmarksId": benchmarks_id,
                        "datasetId": dataset_id,
                        "task_type": task,
                        "tool": method,
                        "Precision-recall AUC": auc_score,
                        "Odds Ratio": oddsratio_score,
                        "created_on": current_date_and_time
                    }
                    if benchmarks_data is not None:
                        labels, y_labels, data = get_bar_plot_data(benchmarks_data, user_results=results)
                        plot_bars(task, labels, y_labels, data)
                    return results
                else: 
                    raise ValueError(f"ccc_pred, ccc_target and score are required for {task}.")

            case "Multimodal Data Integration" | "MI":
                if embedding_key is not None and mod1_key is not None and batch_key is not None and label_key is not None:
                    metrics_dict = multimodal_metrics(mdata, embed=embedding_key, mod1=mod1_key, batch=batch_key, label_key=label_key)
                    results = {**task_info, **metrics_dict}
                    if benchmarks_data is not None:
                        labels, y_labels, data = get_bar_plot_data(benchmarks_data, user_results=results)
                        plot_bars(task, labels, y_labels, data)
                    return results
                else: 
                    raise ValueError(f"embedding_key, mod1_key, label_key and batch_key are required for {task}.")

            case "Cell Type Annotation" | "CT":
                if label_pred_key is not None and label_key is not None:
                    accuracy, f1_macro, f1_micro, f1_weighted = annotation_metrics(adata.obs[label_key], adata.obs[label_pred_key])
                    results = {
                        "benchmarksId": benchmarks_id,
                        "datasetId": dataset_id,
                        "task_type": task,
                        "tool": method,
                        "Accuracy": accuracy,
                        "F1_macro": f1_macro,
                        "F1_micro": f1_micro,
                        "F1_weighted": f1_weighted,
                        "created_on": current_date_and_time
                    }
                    if benchmarks_data is not None:
                        labels, y_labels, data = get_bar_plot_data(benchmarks_data, user_results=results)
                        plot_bars(task, labels, y_labels, data)
                    return results
                else: 
                    raise ValueError(f"label_pred_key, and label_key are required for {task}.")

            case _:  # Default case, equivalent to 'default' in other languages
                raise ValueError(f"{task} is not supported. Please input the task name from the following list [Clustering, Imputation, Batch Integration, Trajectory, Cell-Cell Communication, Multimodal Data Integration, Cell Type Annotation].")
    else:
        raise ValueError("benchmarks_id or task is required.")


def write_json(data, file_path="./output.json"):
    # Open the file in write mode ('w') and use json.dump() to write the dictionary
    with open(file_path, 'w') as json_file:
        json.dump(data, json_file, indent=4, default=serialize_datetime) # indent=4 for pretty-printing

    print(f"Dictionary successfully written to {file_path}")


def serialize_datetime(obj):
    if isinstance(obj, datetime):
        return obj.isoformat()
    raise TypeError("Type not serializable")
