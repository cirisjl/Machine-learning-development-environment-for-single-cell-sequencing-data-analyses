from pymongo import MongoClient
import hashlib

mongo_url = "mongodb://mongodb:65528"
# Connect to MongoDB using the URL
client = MongoClient(mongo_url)
# Access your database and collection
db = client["ai-single-cell"]
datasets_collection = db.datasets
user_datasets_collection = db.get_collection("user-datasets")
pp_results_collection = db.get_collection("pp-results")


def generate_process_id(file_md5, process, method, parameters):
    process_id = None
    if parameters.use_default:
        process_id = hashlib.md5(f"{file_md5}_{process}_{method}_{parameters.assay}".encode("utf_8")).hexdigest()
    else:
        process_id = hashlib.md5(f"{file_md5}_{process}_{method}_{parameters}".encode("utf-8")).hexdigest()
    return process_id


def pp_results_exists(process_id):
    result = pp_results_collection.find_one({'process_id': process_id})
    return result


def create_pp_results(pp_results):
    pp_results_collection.insert_one(pp_results)
    return


def create_datasets(datasets):
    datasets_collection.insert_one(datasets)
    return


def create_user_datasets(datasets):
    user_datasets_collection.insert_one(datasets)
    return