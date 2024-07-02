import hashlib
import os
from pymongo import MongoClient
from pymongo.errors import DuplicateKeyError
from boltons.iterutils import remap


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
    try:
        pp_results_collection.update_one({'process_id': process_id}, {'$set': pp_results}, upsert=True)
    except DuplicateKeyError:
        pp_results_collection.update_one({'process_id': process_id}, {'$set': pp_results})
    if "_id" in pp_results: 
        pp_results.pop("_id")
    return

# Append new process_ids to dataset after each process
def append_pp_ids_to_ds(process_ids, dataset_id):
    collection = datasets_collection
    if dataset_id.split("-")[0] == "U":
        print("Appending new process_ids to user dataset")
        collection = user_datasets_collection
    result = collection.find_one({'Id': dataset_id})
    if  "process_ids" not in result:
        collection.update_one({'Id': dataset_id}, {'$set': {'process_ids': process_ids}})
    else:
        pp_ids = list(set(process_ids)|set(result['process_ids']))
        collection.update_one({'Id': dataset_id}, {'$set': {'process_ids': pp_ids}})
    return


def upsert_jobs(data):
    data = clear_dict(data)
    job_id = data['job_id']
    # data.pop("job_id")
    jobs_collection.update_one({'job_id': job_id}, {'$set': data}, upsert=True)

    if "process_ids" in data and "datasetId" in data: # Append new process_ids to dataset
        append_pp_ids_to_ds(data['process_ids'], data['datasetId'])

    if "_id" in data: 
        data.pop("_id")
    return


def upsert_benchmarks(benchmarksId, results):
    results = clear_dict(results)
    try:
        benchmarks_collection.update_one({'benchmarksId': benchmarksId}, {'$set': results}, upsert=True)
    except DuplicateKeyError:
        benchmarks_collection.update_one({'benchmarksId': benchmarksId}, {'$set': results})
    if "_id" in results: 
        results.pop("_id")
    return


def create_bm_results(process_id, bm_results):
    bm_results = clear_dict(bm_results)
    try:
        bm_results_collection.update_one({'process_id': process_id}, {'$set': bm_results}, upsert=True)
    except DuplicateKeyError:
        bm_results_collection.update_one({'process_id': process_id}, {'$set': bm_results}, upsert=True)
    if "_id" in bm_results: 
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
    if "_id" in results: 
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
