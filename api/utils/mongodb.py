import hashlib
import os
from pymongo import MongoClient
from boltons.iterutils import remap


mongo_url = "mongodb://mongodb:65528"
# Connect to MongoDB using the URL
client = MongoClient(mongo_url)
# Access your database and collection
db = client["ai-single-cell"]
datasets_collection = db.datasets
user_datasets_collection = db.get_collection("user-datasets")
pp_results_collection = db.get_collection("pp-results")
task_results_collection = db.get_collection("task-results")
benchmarks_collection = db.get_collection("benchmarks")
benchmark_results_collection = db.get_collection("benchmark_results")
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


def upsert_task_results(results):
    results = clear_dict(results)
    task_id = results['taskId']
    task_results_collection.update_one({'task_id': task_id}, {'$set': results}, upsert=True)
    if "_id" in results: 
        results.pop("_id")
    return


def upsert_benchmarks(benchmarks_id, results):
    results = clear_dict(results)
    benchmarks_collection.update_one({'benchmarks_id': benchmarks_id}, {'$set': results}, upsert=True)
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


def create_benchmark_results(process_id, benchmark_results):
    benchmark_results = clear_dict(benchmark_results)
    benchmark_results_collection.update_one({'process_id': process_id}, {'$set': benchmark_results}, upsert=True)
    if "_id" in benchmark_results: 
        benchmark_results.pop("_id")
    return


def benchmark_result_exists(process_id):
    result = benchmark_results_collection.find_one({'process_id': process_id}, {'_id': 0})
    return result