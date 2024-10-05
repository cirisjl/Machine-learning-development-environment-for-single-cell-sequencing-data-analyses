from typing import List
from celery import shared_task
from tools.formating.formating import *
from tools.run_qc import run_qc
from tools.run_normalization import run_normalization
from tools.run_imputation import run_imputation
from tools.run_integration import run_integration
from tools.run_evaluation import run_evaluation
from tools.run_reduction import run_reduction
from tools.run_conversion import run_conversion
from tools.load_metadata import load_metadata
from benchmarks.run_benchmarks import run_benchmarks
from benchmarks.run_data_split import run_data_split
from benchmarks.run_subset_data import run_subset_data
from workflows.clustering import run_clustering


@shared_task(bind=True, name='tools:create_qc_task') 
def create_qc_task(self, ds_dict:dict):
    job_id = self.request.id
    results = run_qc(job_id, ds_dict)
    return results


@shared_task(bind=True, name='tools:create_normalization_task') 
def create_normalization_task(self, ds_dict:dict):
    job_id = self.request.id
    results = run_normalization(job_id, ds_dict)
    return results


@shared_task(bind=True, name='tools:create_imputation_task') 
def create_imputation_task(self, ds_dict:dict):
    job_id = self.request.id
    results = run_imputation(job_id, ds_dict)
    return results


@shared_task(bind=True, name='tools:create_reduction_task') 
def create_reduction_task(self, ds_dict:dict):
    job_id = self.request.id
    results = run_reduction(job_id, ds_dict)
    return results


@shared_task(bind=True, name='tools:create_conversion_task') 
def create_conversion_task(self, ds_dict:dict):
    job_id = self.request.id
    results = run_conversion(job_id, ds_dict)
    return results


@shared_task(bind=True, name='tools:create_integration_task') 
def create_integration_task(self, ids_dict:dict):
    job_id = self.request.id
    results = run_integration(job_id, ids_dict)
    return results


@shared_task(bind=True, name='tools:load_metadata_task') 
def load_metadata_task(self, file_dict:dict):
    job_id = self.request.id
    results = load_metadata(job_id, file_dict)
    return results


@shared_task(bind=True, name='tools:create_evaluation_task') 
def create_evaluation_task(self, dataset, input, userID, output, methods, layer=None, genes=None, ncores=12, show_error=True):
    job_id = self.request.id
    results = run_evaluation(job_id, dataset, input, userID, output, methods, layer=None, genes=None, ncores=12, show_error=True)
    return results


# Benchmarks
@shared_task(bind=True, name='benchmarks:create_benchmarks_task') 
def create_benchmarks_task(self, task_dict:dict):
    job_id = self.request.id
    results = run_benchmarks(job_id, task_dict)
    return results


@shared_task(bind=True, name='benchmarks:create_data_split_task') 
def create_data_split_task(self, task_dict:dict):
    job_id = self.request.id
    results = run_data_split(job_id, task_dict)
    return results


@shared_task(bind=True, name='benchmarks:create_subset_data_task') 
def create_subset_data_task(self, task_dict:dict):
    job_id = self.request.id
    results = run_subset_data(job_id, task_dict)
    return results


# Workflows
@shared_task(bind=True, name='workflows:create_clustering_task') 
def create_clustering_task(self, ds_dict:dict):
    job_id = self.request.id
    results = run_clustering(job_id, ds_dict)
    return results