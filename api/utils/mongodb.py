import hashlib
import os
from pymongo import MongoClient
from boltons.iterutils import remap


mongo_url = "mongodb://mongodb:65528"
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
    pp_results_collection.update_one({'process_id': process_id}, {'$set': pp_results}, upsert=True)
    if "_id" in pp_results: 
        pp_results.pop("_id")
    return


def upsert_jobs(updates):
    updates = clear_dict(updates)
    job_id = updates['job_id']
    jobs_collection.update_one({'job_id': job_id}, {'$set': updates}, upsert=True)
    if "_id" in updates: 
        updates.pop("_id")
    return


def upsert_benchmarks(benchmarksId, results):
    results = clear_dict(results)
    benchmarks_collection.update_one({'benchmarksId': benchmarksId}, {'$set': results}, upsert=True)
    if "_id" in results: 
        results.pop("_id")
    return


def upsert_workflows(workflows_id, results):
    results = clear_dict(results)
    workflows_collection.update_one({'workflows_id': workflows_id}, {'$set': results}, upsert=True)
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


def create_bm_results(process_id, bm_results):
    bm_results = clear_dict(bm_results)
    bm_results_collection.update_one({'process_id': process_id}, {'$set': bm_results}, upsert=True)
    if "_id" in bm_results: 
        bm_results.pop("_id")
    return


def benchmark_result_exists(process_id):
    result = bm_results_collection.find_one({'process_id': process_id}, {'_id': 0})
    return result