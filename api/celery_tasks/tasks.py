from typing import List
from celery import shared_task
from tools.run_qc import run_qc
from tools.run_normalization import run_normalization
from tools.run_imputation import run_imputation

@shared_task(bind=True,autoretry_for=(Exception,), retry_backoff=True, retry_kwargs={"max_retries": 5},
             name='tools:create_qc_task')
def create_qc_task(self, dataset, input, userID, output, methods, path_of_scrublet_calls='./scrublet_calls.tsv', show_error=True):
    results = run_qc(dataset, input,userID, output, methods, show_error = True)
    return results


@shared_task(bind=True,autoretry_for=(Exception,), retry_backoff=True, retry_kwargs={"max_retries": 5},
             name='tools:create_normalization_task')
def create_normalization_task(self, dataset, input, userID, output, methods, default_assay='RNA', output_format='AnnData', species=None, idtype='ENSEMBL',show_umap = True, show_error = True):
    results = run_normalization(dataset, input, userID, output, methods, default_assay='RNA', output_format='AnnData', species=None, idtype='ENSEMBL', show_error = True)
    return results


@shared_task(bind=True, autoretry_for=(Exception,), retry_backoff=True, retry_kwargs={"max_retries": 5},
             name='tools:create_imputation_task')
def create_imputation_task(self, dataset, input, userID, output, methods, layer=None, genes=None, ncores=12, show_error=True):
    results = run_imputation(dataset, input, userID,  output, methods, layer=None, genes=None, ncores=12, show_error=True)
    return results
