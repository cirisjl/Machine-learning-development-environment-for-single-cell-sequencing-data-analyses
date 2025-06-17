import hashlib
import os
from pymongo import MongoClient
from pymongo.errors import DuplicateKeyError
from boltons.iterutils import remap
from tools.utils.gzip_str import *
from tools.formating.formating import regularise_df
import json_numpy


mongo_url = "mongodb://mongodb:65530"
# Connect to MongoDB using the URL
client = MongoClient(mongo_url)
# Access your database and collection
db = client["oscb"]
datasets_collection = db.datasets
user_datasets_collection = db.get_collection("user_datasets")
pp_results_collection = db.get_collection("pp_results")
jobs_collection = db.get_collection("jobs")
benchmarks_collection = db.get_collection("benchmarks")
bm_results_collection = db.get_collection("bm_results")
workflows_collection = db.get_collection("workflows")

pp_results_collection.create_index({'process_id': 1}, unique=True, background=True)
bm_results_collection.create_index({'process_id': 1}, unique=True, background=True)
workflows_collection.create_index({'workflows_id': 1}, unique=True, background=True)
datasets_collection.create_index({'Id': 1}, unique=True, background=True)
benchmarks_collection.create_index({'benchmarksId': 1}, unique=True, background=True)
jobs_collection.create_index({'job_id': 1}, unique=True, background=True)

def generate_process_id(file_md5, process, method, parameters):
    if isinstance(method, list) and len(method) == 1: method = method[0] # If there is only 1 item in a list, then it should be taken as a string.
    process_id = hashlib.md5(f"{file_md5}_{process}_{method}_{parameters}".encode("utf-8")).hexdigest()
    return process_id


def generate_workflow_id(file_md5, workflow, parameters):
    workflow_id = hashlib.md5(f"{file_md5}_{workflow}_{parameters}".encode("utf-8")).hexdigest()
    return workflow_id


def pp_result_exists(process_id):
    result = pp_results_collection.find_one({'process_id': process_id}, {'_id': 0})
    return result


def create_pp_results(process_id, pp_results):
    pp_results = clear_dict(pp_results)
    # pp_results = removeNullNoneEmpty(pp_results)
    if "obs" in pp_results.keys():
        pp_results.pop("obs")  # Remove obs key if it exists, as it is not needed in the database

    try:
        pp_results_collection.update_one({'process_id': process_id}, {'$set': pp_results}, upsert=True)
    except DuplicateKeyError:
        pp_results_collection.update_one({'process_id': process_id}, {'$set': pp_results})
    if "_id" in pp_results.keys(): 
        pp_results.pop("_id")
    return


def get_pp_results(process_ids, umap=False, record_type=None):
    pp_results = None
    if not umap:
        if record_type == None:
            pp_results = pp_results_collection.find({'process_id': { "$in": process_ids }}, { "_id": 0, "process_id": 1, "description": 1, "parameters": 1, "stage": 1, "process": 1, "method": 1, "nCells": 1, "adata_path": 1, "md5": 1, "info": 1, "cell_metadata": 1, "obs_names": 1, "default_assay": 1, "assay_names": 1, "umap": 1, "umap_3d": 1, "highest_expr_genes": 1, "evaluation_results": 1 ,"layers" : 1, "layer" : 1, "tsne": 1, "tsne_3d": 1})
        elif record_type == 'table':
            pp_results = pp_results_collection.find({'process_id': { "$in": process_ids }}, { "_id": 0, "process_id": 1, "description": 1, "stage": 1, "process": 1, "method": 1, "nCells": 1, "adata_path": 1, "md5": 1, "info": 1, "cell_metadata": 1, "obs_names": 1, "default_assay": 1, "assay_names": 1,"layers" : 1, "layer" : 1})
    else:
        pp_results = pp_results_collection.find({'process_id': {"$in": process_ids }}, { "_id": 0, "tsne": 0, "process_id": 0, "description": 0, "stage": 0, "process": 0, "method": 0, "nCells": 0, "adata_path": 0, "md5": 0, "info": 0, "default_assay": 0, "assay_names": 0, "highest_expr_genes": 0, "evaluation_results": 0 })
    pp_results = list(pp_results)
    results = []

    if len(pp_results) != 0:
        for pp_result in pp_results:
            if 'cell_metadata' in pp_result.keys():
                pp_result['obs'] = pp_result['cell_metadata']
                obs_dict = gunzip_dict(pp_result['cell_metadata'])
                obs = pd.DataFrame.from_dict(obs_dict)
                obs = obs.set_index('index')  
                pp_result['cell_metadata'] = obs

            if 'genes' in pp_result.keys():
                genes = gunzip_list(pp_result['genes'])

            # if 'gene_metadata' in pp_result.keys():
            #     pp_result['gene_metadata'] = gunzip_df(pp_result['gene_metadata'])

            if 'umap' in pp_result.keys():
                pp_result['umap'] = json_numpy.loads(pp_result['umap'])

            if 'umap_3d' in pp_result.keys():
                pp_result['umap_3d'] =json_numpy.loads(pp_result['umap_3d'])

            if 'tsne' in pp_result.keys():
                pp_result['tsne'] = json_numpy.loads(pp_result['tsne'])
            
            if 'highest_expr_genes' in pp_result.keys():
                pp_result['highest_expr_genes']['counts_top_genes'] = json_numpy.loads(pp_result['highest_expr_genes']['counts_top_genes'])
            
            results.append(pp_result)
    
    return results


# Append new process_ids to dataset after each process
def append_pp_ids_to_ds(process_ids, datasetId):
    collection = datasets_collection
    if datasetId.split("-")[0] == "U":
        print("Appending new process_ids to user dataset")
        collection = user_datasets_collection
    result = collection.find_one({'Id': datasetId})
    if  "process_ids" not in result.keys():
        collection.update_one({'Id': datasetId}, {'$set': {'process_ids': process_ids}})
    else:
        pp_ids = list(set(process_ids)|set(result['process_ids']))
        collection.update_one({'Id': datasetId}, {'$set': {'process_ids': pp_ids}})
    return


def upsert_jobs(data):
    data = clear_dict(data)
    job_id = data['job_id']
    # data.pop("job_id")
    jobs_collection.update_one({'job_id': job_id}, {'$set': data}, upsert=True)

    if "process_ids" in data.keys():
        if "datasetId" in data.keys(): # Append new process_ids to dataset
            append_pp_ids_to_ds(data['process_ids'], data['datasetId'])
        if "datasetIds" in data.keys(): # Append new process_ids to dataset
            for datasetId in data['datasetIds']:
                append_pp_ids_to_ds(data['process_ids'], datasetId)

    if "_id" in data.keys(): 
        data.pop("_id")
    return


def upsert_benchmarks(benchmarksId, results):
    results = clear_dict(results)
    try:
        benchmarks_collection.update_one({'benchmarksId': benchmarksId}, {'$set': results}, upsert=True)
    except DuplicateKeyError:
        benchmarks_collection.update_one({'benchmarksId': benchmarksId}, {'$set': results})
    if "_id" in results.keys(): 
        results.pop("_id")
    return


def create_bm_results(process_id, bm_results):
    bm_results = clear_dict(bm_results)
    try:
        bm_results_collection.update_one({'process_id': process_id}, {'$set': bm_results}, upsert=True)
    except DuplicateKeyError:
        bm_results_collection.update_one({'process_id': process_id}, {'$set': bm_results}, upsert=True)
    if "_id" in bm_results.keys(): 
        bm_results.pop("_id")
    return


def benchmark_result_exists(process_id):
    result = bm_results_collection.find_one({'process_id': process_id}, {'_id': 0})
    return result


def upsert_workflows(workflows_id, results):
    results = clear_dict(results)
    try:
        workflows_collection.update_one({'workflows_id': workflows_id}, {'$set': results}, upsert=True)
    except DuplicateKeyError:
        workflows_collection.update_one({'workflows_id': workflows_id}, {'$set': results})
    if "_id" in results.keys(): 
        results.pop("_id")
    return


def create_datasets(datasets):
    datasets_collection.insert_one(datasets)
    return


def create_user_datasets(datasets):
    user_datasets_collection.insert_one(datasets)
    return


def clear_dict(d):
    drop_falsey = lambda path, key, value: value is not None and value != [] and value != {} and value != [{}]
    d = remap(d, visit=drop_falsey)
    return d


def removeNullNoneEmpty(ob):
    l = {}
    for k, v in ob.items():
        if(isinstance(v, dict)):
            x = removeNullNoneEmpty(v)
            if(len(x.keys())>0):
                l[k] = x
        
        elif(isinstance(v, list)):
            p = []
            for c in v:
                if(isinstance(c, dict)):
                    x = removeNullNoneEmpty(c)
                    if(len(x.keys())>0):
                        p.append(x)
                elif(c is not None and c != ''):
                    p.append(c)
            l[k] = p
        elif(v is not None and v!=''):
            l[k] = v
    return l


