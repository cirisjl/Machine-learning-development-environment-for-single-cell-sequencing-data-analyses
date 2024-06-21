from tools.formating.formating import *
from tools.run_qc import run_qc
from tools.run_normalization import run_normalization
from tools.run_imputation import run_imputation
from utils.unzip import unzip_file_if_compressed
from fastapi import HTTPException, status
from utils.redislogger import *
from utils.mongodb import generate_workflow_id, upsert_jobs, upsert_workflows
from datetime import datetime

def run_clustering(job_id, ds:dict, random_state=0):
    wf_results = {}
    process_ids = []
    userID = ds['userID']
    output = None
    input = unzip_file_if_compressed(job_id, ds['input'])
    md5 = get_md5(input)
    workflow_id = generate_workflow_id(md5, "clustering", ds)

    # wf_results = pp_result_exists(process_id)
    upsert_jobs(
        {
            "job_id": job_id, 
            "created_by": userID,
            "completed_on": datetime.now(),
            "status": "Processing"
        }
    )

    qc_results = run_qc(job_id, ds)
    if qc_results is not None:
        ds['input'] = qc_results['output']
        wf_results['QC'] = qc_results["process_ids"]
        wf_results['QC_ouput'] = qc_results['ouput']
        if ds["normalization_params"]["methods"] is not None:
            normalization_results = run_normalization(job_id, ds)
            wf_results['normalization'] = normalization_results["process_ids"]
            wf_results['normalization_ouput'] = normalization_results['ouput']
            output = normalization_results['ouput']
        elif ds["imputationParameters"]["methods"] is not None:
            imputation_results = run_imputation(job_id, ds)
            wf_results['imputation'] = imputation_results["process_ids"]
            wf_results['imputation_ouput'] = imputation_results['ouput']
            output = imputation_results['ouput']
        upsert_workflows(workflow_id, wf_results)

    results = {
        "output": output,
        "workflow_id": workflow_id
    }
    
    upsert_jobs(
        {
            "job_id": job_id, 
            "output": output,
            "workflow_id": workflow_id,
            "results": results,
            "completed_on": datetime.now(),
            "status": "Success"
        }
    )

    return results

        


    