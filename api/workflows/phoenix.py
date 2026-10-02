import os
import subprocess
from fastapi import HTTPException, status
from tools.formating.formating import get_md5
from utils.redislogger import *
from utils.mongodb import generate_workflow_id, upsert_jobs, upsert_workflows
from utils.wgs import create_samplesheet
from datetime import datetime
from exceptions.custom_exceptions import CeleryTaskException


def run_phoenix_wf(job_id, ds:dict, random_state=0):
    wf_results = {}
    mode = ds['mode']
    wgs_input = ds['input']
    kraken2db = ds['kraken2db']
    scaffold_ext = ds['scaffold_ext']
    input_sra = ds['input_sra']
    use_sra = ds['use_sra']
    centar = ds['centar']
    griphin_out = ds['griphin_out']
    outdir = ds['outdir']
    coverage = ds['coverage']
    busco_db_path = ds['busco_db_path']
    create_ncbi_sheet = ds['create_ncbi_sheet']
    email = ds['email']
    multiqc_title = ds['multiqc_title']

    description = ds['description']
    process_ids = []
    userID = ds['userID']
    fig_path = None
    method = ds['method']
    process = ds['process']
    max_cpus = os.environ.get("MAX_CPUS", 16)
    max_memory = os.environ.get("MAX_MEMORY", '128.GB')
    max_time = os.environ.get("MAX_TIME", '240.h')
    phoenix_version = os.environ.get("PHOENIX_VERSION", "v2.3.2")
    kraken2db_home = os.environ.get("KRAKEN2DB_HOME", "/ps/phoenix/kraken2db/")
    md5 = get_md5(wgs_input)

    kraken2db = os.path.join(kraken2db_home, kraken2db)
    if not os.path.exists(kraken2db):
        raise CeleryTaskException(f"Kraken2DB does not exist!")

    if not os.path.exists(wgs_input):
        raise CeleryTaskException(f"File or folder does not exist!")
    if os.path.isdir(wgs_input):
        wgs_input = create_samplesheet(wgs_input)
        if outdir is None:
            outdir = os.path.join(wgs_input, 'phoenix')
    elif not os.path.isfile(wgs_input) or not wgs_input.endswith(".csv"):
        raise CeleryTaskException(f"Input is not a valid CSV file!")
    elif outdir is None:
        outdir = os.path.join(os.path.dirname(wgs_input), 'phoenix')

    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)

    if griphin_out is not None and not os.path.exists(griphin_out):
        os.makedirs(griphin_out, exist_ok=True)

    # Create folders for output figures
    workflow_id = generate_workflow_id(md5, "PHoeNIx", ds)
    fig_path = os.path.join(outdir, 'figures')
    if not os.path.exists(fig_path):
        os.makedirs(fig_path, exist_ok=True)
    wf_results['figures'] = fig_path

    upsert_jobs(
        {
            "job_id": job_id, 
            "created_by": userID,
            "workflow_id": workflow_id,
            "Description": description,
            "Method": method,
            "datasetURL": input,
            "datasetId": md5,
            "Process": "PHoeNIx",
            "Category": 'WGS Workflow',
            "Created on": datetime.now(),
            "Status": "Processing"
        }
    )

    try:
        report_path = os.path.join(outdir, 'Phoenix_Summary.tsv')
        cmd = f"nextflow run cdcgov/phoenix -r {phoenix_version} -profile docker --mode {mode} --kraken2db {kraken2db} --input {wgs_input} --outdir {outdir} --coverage {coverage} --max_cpus {max_cpus} --max_memory {max_memory} --max_time {max_time}"

        if (mode == "SRA" or mode == "CDC_SRA") and use_sra and input_sra is not None:
            cmd += f" --input_sra {input_sra} --use_sra"

        if mode == "CENTAR" and centar:
            cmd += f" --centar"

        if create_ncbi_sheet:
            cmd += f" --create_ncbi_sheet"
        
        if email is not None:
            cmd += f" --email {email}"

        # if griphin_out is not None:
        #     cmd += f" --griphin_out {griphin_out}"

        # if mode == "combine_griphins" and griphin_out is not None:
        #     cmd += f" --griphin_out {griphin_out}"

        if multiqc_title is not None:
            cmd += f" --multiqc_title {multiqc_title}"

        if (mode == "SCAFFOLDS" or mode == "CDC_SCAFFOLDS") and scaffolds_ext is not None:
            cmd += f" --scaffolds_ext {scaffolds_ext}"

        # Execute the Nextflow command
        s = subprocess.call(cmd, shell=True)

        if os.path.exists(report_path):
            results = {
                "output": outdir,
                "workflow_id": workflow_id,
                "md5": md5,
                "wf_results": report_path,
                # "figures":fig_path, 
                "workflow_id": workflow_id,
                # "process_ids": process_ids
            }
            
            upsert_jobs(
                {
                    "job_id": job_id, 
                    "output": outdir,
                    # "process_ids": process_ids,
                    "workflow_id": workflow_id,
                    "results": results,
                    # "figures": fig_path, 
                    "Completed on": datetime.now(),
                    "Status": "Success"
                }
            )
        else:
            detail = f'PHoeNIx workflow is failed: Report not generated.'
            redislogger.error(job_id, detail)
            upsert_jobs(
                {
                    "job_id": job_id, 
                    "results": {"error": detail},
                    "Completed on": datetime.now(),
                    "Status": "Failure"
                }
            )
            raise CeleryTaskException(detail)
    except Exception as e:
        detail = f"PHoeNIx workflow is failed: {e}"
        upsert_jobs(
            {
                "job_id": job_id, 
                "results": {"error": detail},
                "Completed on": datetime.now(),
                "Status": "Failure"
            }
        )
        raise CeleryTaskException(detail)

    return results 