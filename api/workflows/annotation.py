import os
from tools.formating.formating import *
from tools.run_qc import run_qc
from tools.run_integration import run_integration
from tools.run_annotation import run_annotation
from utils.unzip import unzip_file_if_compressed
from fastapi import HTTPException, status
from utils.redislogger import *
from utils.mongodb import generate_workflow_id, upsert_jobs, upsert_workflows
from datetime import datetime

def run_annotation_wf(job_id, dss:dict, random_state=0):
    wf_results = {}
    process_ids = []
    userID = dss['userID']
    datasetIds = dss['datasetIds']
    datasets = dss['dataset']
    inputs = dss['input']
    do_umap = dss['do_umap']
    do_cluster = dss['do_cluster']
    user_refs = dss['user_refs']
    qc_params = dss['qc_params']
    integration_params = dss['integration_params']
    annotation_params = dss['annotation_params']
    reduction_params = dss['reduction_params']
    n_neighbors = reduction_params['n_neighbors']
    n_pcs = reduction_params['n_pcs']
    resolution = reduction_params['resolution']
    qc_params['n_neighbors'] = n_neighbors
    qc_params['n_pcs'] = n_pcs
    qc_params['resolution'] = resolution
    integration_params['dims'] = n_neighbors
    integration_params['npcs'] = n_pcs 
    integration_params['resolution'] = resolution
    annotation_params['n_neighbors'] = n_neighbors
    annotation_params['n_pcs'] = n_pcs 
    annotation_params['resolution'] = resolution
    fig_path = None
    md5 = []
    adata_outputs = None
    final_outputs = []

    # Initialize methodMap
    methodMap = {}

    # Extract methods for Quality Control
    if qc_params and "methods" in qc_params and qc_params["methods"]:
        methodMap["Quality Control"] = ", ".join(qc_params["methods"])

    # Extract methods for Integration
    if integration_params and "methods" in integration_params and integration_params["methods"]:
        methodMap["Integration"] = ", ".join(integration_params["methods"])

    # Extract methods for Annotation
    if annotation_params and "methods" in annotation_params and annotation_params["methods"]:
        methodMap["Annotation"] = ", ".join(annotation_params["methods"])
        
    output = None
    description = "Annotation Workflow"

    abs_inputList = []

    if len(inputs) > 0:
        for input in inputs:
            if input is not None:
                input = unzip_file_if_compressed(job_id, input)
                md5 = md5 + get_md5(input)
                abs_inputList.append(input)
    else:
        raise CeleryTaskException("No input file is found.")

    # Create folders for output figures
    workflow_id = generate_workflow_id(md5, "annotation", dss)
    fig_path = os.path.join(os.path.dirname(inputs[0]), 'workflow', workflow_id)
    if not os.path.exists(fig_path):
        os.makedirs(fig_path, exist_ok=True)
    wf_results['figures'] = fig_path

    # for method i want methodMap key value pairs.

    if datasetIds is not None:
        description = f"Annotation Workflow for {datasetIds}"
    elif dataset is not None:
        description = f"Annotation Workflow for {datasets}"

    # wf_results = pp_result_exists(process_id)
    upsert_jobs(
        {
            "job_id": job_id, 
            "created_by": userID,
            "description": description,
            "method": str(methodMap).replace("'", "").replace("{", "").replace("}", ""),
            "datasetURL": str(inputs).replace("'", "").replace("{", "").replace("}", ""),
            "datasetId": str(datasetIds).replace("'", "").replace("{", "").replace("}", ""),
            "process": "Annotation",
            "category": 'workflow',
            "created_on": datetime.now(),
            "status": "Processing"
        }
    )

    integration_inputs = []
    qc_process_ids = []
    for i in range(len(abs_inputList)):
        ds = {}
        ds['userID'] = userID
        ds['input'] = abs_inputList[i]
        ds['output'] = dss['output']
        ds['datasetId'] = datasetIds[i]
        ds['dataset'] = datasets[i]
        ds['do_umap'] = do_umap
        ds['do_cluster'] = do_cluster
        ds['qc_params'] = qc_params
        
        qc_results = run_qc(job_id, ds, fig_path=fig_path)
        if qc_results is not None:
            integration_inputs.append(qc_results['adata_path'])
            process_ids.extend(qc_results["process_ids"])
            qc_process_ids.extend(qc_results["process_ids"])

    wf_results['QC'] = qc_process_ids
    wf_results['QC_output'] = integration_inputs

    if len(integration_inputs) > 0 and len(integration_params["methods"]) > 0:
        dss['input'] = integration_inputs
        integration_results = run_integration(job_id, dss, fig_path=fig_path)
        wf_results['integration'] = integration_results["process_ids"]
        process_ids.extend(integration_results["process_ids"])
        wf_results['integration_output'] = integration_results['output']
        output = integration_results['output']
        adata_outputs = integration_results['adata_path']
    else:
        raise CeleryTaskException("No integration input file is found.")
    
    if len(adata_outputs) > 0 and len(annotation_params["methods"]) > 0:
        ann_process_ids = []
        ann_outputs = []
        for adata_output in adata_outputs:
            ds = {}
            ds['userID'] = userID
            ds['input'] = adata_output
            ds['output'] = dss['output']
            ds['datasetId'] = datasetIds[0]
            ds['dataset'] = datasets[0]
            ds['species'] = dss['species']
            ds['user_refs'] = user_refs
            ds['do_umap'] = do_umap
            ds['do_cluster'] = do_cluster
            ds['annotation_params'] = annotation_params

            annotation_results = run_annotation(job_id, ds, fig_path=fig_path)
            ann_process_ids.extend(annotation_results["process_ids"])
            process_ids.extend(annotation_results["process_ids"])
            ann_outputs = ann_outputs + annotation_results['output']
            final_outputs.append(annotation_results['adata_path'])

        wf_results['annotation'] = ann_process_ids
        wf_results['annotation_output'] = ann_outputs
        output = ann_outputs
    else:
        raise CeleryTaskException("No annotation input file is found.")

    results = {
        "output": output,
        # "workflow_id": workflow_id,
        "md5": md5,
        "adata_path": final_outputs,
        "wf_results": wf_results,
        # "figures":fig_path, 
        "process_ids": process_ids
    }
    
    upsert_jobs(
        {
            "job_id": job_id, 
            "output": output,
            "process_ids": process_ids,
            # "workflow_id": workflow_id,
            "results": results,
            # "figures": fig_path, 
            "completed_on": datetime.now(),
            "status": "Success"
        }
    )

    return results