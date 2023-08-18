from typing import List
from celery import shared_task
from tools.run_qc import run_qc
from tools.run_normalization import run_normalization
from tools.run_imputation import run_imputation
from tools.run_integration import run_integration


@shared_task(bind=True,autoretry_for=(Exception,), retry_backoff=True, retry_kwargs={"max_retries": 5},
             name='tools:create_qc_task')
def create_qc_task(self, dataset, input, userID, output, methods, path_of_scrublet_calls='./scrublet_calls.tsv', show_error=True):
    task_id = self.request.id
    results = run_qc(task_id,dataset, input,userID, output, methods, show_error = True)
    return results


@shared_task(bind=True,autoretry_for=(Exception,), retry_backoff=True, retry_kwargs={"max_retries": 5},
             name='tools:create_normalization_task')
def create_normalization_task(self, dataset, input, userID, output, methods,species, default_assay='RNA', output_format='AnnData', idtype='ENSEMBL',show_umap = True, show_error = True):
    task_id = self.request.id
    results = run_normalization(task_id, dataset, input, userID, output, methods,species, default_assay='RNA', output_format='AnnData', idtype='ENSEMBL', show_error = True)
    return results


@shared_task(bind=True, autoretry_for=(Exception,), retry_backoff=True, retry_kwargs={"max_retries": 5},
             name='tools:create_imputation_task')
def create_imputation_task(self, dataset, input, userID, output, methods, layer=None, genes=None, ncores=12, show_error=True):
    task_id = self.request.id
    results = run_imputation(task_id,dataset, input, userID,  output, methods, layer=None, genes=None, ncores=12, show_error=True)
    return results

@shared_task(bind=True,autoretry_for=(Exception,), retry_backoff=True, retry_kwargs={"max_retries": 5},
             name='tools:create_integration_task')
def create_integration_task(self, datasets, inputs, userID, output, methods,species, default_assay='RNA', output_format='Seurat', genes=None, reference=12, show_error = True):
    task_id = self.request.id
    results = run_integration(task_id, datasets, inputs, userID, output, methods, species, default_assay='RNA', output_format='Seurat', genes=None, reference=12, show_error = True)
    return results

@shared_task(bind=True, autoretry_for=(Exception,), retry_backoff=True, retry_kwargs={"max_retries": 5},
             name='tools:create_evaluation_task')
def create_evaluation_task(self, dataset, input, userID, output, methods, layer=None, genes=None, ncores=12, show_error=True):
    task_id = self.request.id
    results = run_evaluation(task_id,dataset, input, userID,  output, methods, layer=None, genes=None, ncores=12, show_error=True)
    return results
