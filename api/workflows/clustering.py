from tools.formating.formating import *
from tools.run_qc import run_qc
from tools.run_normalization import run_normalization
from tools.run_imputation import run_imputation
from utils.unzip import unzip_file_if_compressed
from fastapi import HTTPException, status
from utils.redislogger import *
from utils.mongodb import generate_workflow_id, upsert_task_results, upsert_workflows

def run_clustering(task_id, ds:dict, random_state=0):
    wf_results = {}
    process_ids = []
    userID = ds['userID']
    output = None
    input = unzip_file_if_compressed(task_id, ds['input'])
    md5 = get_md5(input)
    workflow_id = generate_workflow_id(md5, "clustering", ds)

    # wf_results = pp_result_exists(process_id)

    qc_results = run_qc(task_id, ds)
    if qc_results is not None:
        ds['input'] = qc_results['output']
        wf_results['QC'] = qc_results["process_ids"]
        wf_results['QC_ouput'] = qc_results['ouput']
        if ds["normalization_params"]["methods"] is not None:
            normalization_results = run_normalization(task_id, ds)
            wf_results['normalization'] = normalization_results["process_ids"]
            wf_results['normalization_ouput'] = normalization_results['ouput']
            output = normalization_results['ouput']
        elif ds["imputationParameters"]["methods"] is not None:
            imputation_results = run_imputation(task_id, ds)
            wf_results['imputation'] = imputation_results["process_ids"]
            wf_results['imputation_ouput'] = imputation_results['ouput']
            output = imputation_results['ouput']
        upsert_workflows(workflow_id, wf_results)

    results = {
        "taskId": task_id,
        "owner": userID,
        "inputfile": input,
        "output": output,
        "workflow_id": workflow_id,
        "status": "Success"
    }
    
    upsert_task_results(results)

    return results

        


    