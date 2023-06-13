from typing import List
from celery import shared_task
from api import tools

@shared_task(bind=True,autoretry_for=(Exception,), retry_backoff=True, retry_kwargs={"max_retries": 5},
             name='tools:create_qc_task')
def create_qc_task(self, dataset, input, userID, output, methods, path_of_scrublet_calls='./scrublet_calls.tsv', show_error=True):
    results = tools.qc(dataset, input,userID, output, methods, path_of_scrublet_calls='./scrublet_calls.tsv', show_error = True)
    return results


@shared_task(bind=True,autoretry_for=(Exception,), retry_backoff=True, retry_kwargs={"max_retries": 5},
             name='tools:create_normalization_task')
def create_normalization_task(self, dataset, input, output, methods, default_assay='RNA', output_format='AnnData', species=None, idtype='ENSEMBL', show_error = True):
    results = tools.normalize(dataset, input, output, methods, default_assay='RNA', output_format='AnnData', species=None, idtype='ENSEMBL', show_error = True)
    return results


@shared_task(bind=True, autoretry_for=(Exception,), retry_backoff=True, retry_kwargs={"max_retries": 5},
             name='tools:create_imputation_task')
def create_imputation_task(self, dataset, input, output, methods, layer=None, genes=None, ncores=12, show_error=True):
    results = tools.impute(dataset, input, output, methods, layer=None, genes=None, ncores=12, show_error=True)
    return results
