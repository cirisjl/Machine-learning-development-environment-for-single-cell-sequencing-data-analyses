import os
import subprocess
from fastapi import HTTPException, status
from tools.formating.formating import get_md5
from utils.redislogger import *
from utils.mongodb import generate_workflow_id, upsert_jobs, upsert_workflows
from utils.wgs import create_samplesheet
from datetime import datetime
from exceptions.custom_exceptions import CeleryTaskException


def run_mycosnp_wf(job_id, ds:dict, random_state=0):
    wf_results = {}
    workflow = ds['workflow']
    wgs_input = ds['input']
    add_sra_file = ds['add_sra_file']
    add_vcf_file = ds['add_vcf_file']
    outdir = ds['outdir']

    assembler = ds['pre_mycosnp_params']['assembler']
    shovill_depth = ds['pre_mycosnp_params']['shovill_depth']
    genome_size = ds['pre_mycosnp_params']['genome_size']
    min_contig_cov = ds['pre_mycosnp_params']['min_contig_cov']
    min_contig_len = ds['pre_mycosnp_params']['min_contig_len']
    gambit_db = ds['pre_mycosnp_params']['gambit_db']
    gambit_h5_dir = ds['pre_mycosnp_params']['gambit_h5_dir']
    subtype_db = ds['pre_mycosnp_params']['subtype_db']

    fasta = ds['nfcore_mycosnp_params']['fasta']
    ref_dir = ds['nfcore_mycosnp_params']['ref_dir']
    ref_masked_fasta = ds['nfcore_mycosnp_params']['ref_masked_fasta']
    ref_fai = ds['nfcore_mycosnp_params']['ref_fai']
    ref_dict = ds['nfcore_mycosnp_params']['ref_dict']
    ref_bwa = ds['nfcore_mycosnp_params']['ref_bwa']
    snpeffconfig = ds['nfcore_mycosnp_params']['snpeffconfig']
    species = ds['nfcore_mycosnp_params']['species']
    positions = ds['nfcore_mycosnp_params']['positions']
    genes = ds['nfcore_mycosnp_params']['genes']
    exclude = ds['nfcore_mycosnp_params']['exclude']
    save_reference = ds['nfcore_mycosnp_params']['save_reference']
    save_alignment = ds['nfcore_mycosnp_params']['save_alignment']
    rapidnj = ds['nfcore_mycosnp_params']['rapidnj']
    fasttree = ds['nfcore_mycosnp_params']['fasttree']
    snpeff = ds['nfcore_mycosnp_params']['snpeff']
    iqtree = ds['nfcore_mycosnp_params']['iqtree']
    raxmlng = ds['nfcore_mycosnp_params']['raxmlng']
    save_debug = ds['nfcore_mycosnp_params']['save_debug']
    skip_samples = ds['nfcore_mycosnp_params']['skip_samples']
    skip_samples_file = ds['nfcore_mycosnp_params']['skip_samples_file']
    skip_combined_analysis = ds['nfcore_mycosnp_params']['skip_combined_analysis']
    skip_phylogeny = ds['nfcore_mycosnp_params']['skip_phylogeny']
    sample_ploidy = ds['nfcore_mycosnp_params']['sample_ploidy']
    coverage = ds['nfcore_mycosnp_params']['coverage']
    gvcfs_filter = ds['nfcore_mycosnp_params']['gvcfs_filter']
    gatkgenotypes_filter = ds['nfcore_mycosnp_params']['gatkgenotypes_filter']
    max_amb_samples = ds['nfcore_mycosnp_params']['max_amb_samples']
    max_perc_amb_samples = ds['nfcore_mycosnp_params']['max_perc_amb_samples']
    min_depth = ds['nfcore_mycosnp_params']['min_depth']
    mask = ds['nfcore_mycosnp_params']['mask']
    
    description = ds['description']
    # process_ids = []
    userID = ds['userID']
    fig_path = None
    method = ds['method']
    process = ds['process']
    max_cpus = os.environ.get("MAX_CPUS", 16)
    max_memory = os.environ.get("MAX_MEMORY", '128.GB')
    max_time = os.environ.get("MAX_TIME", '240.h')
    mycosnp_version = os.environ.get("MYCOSNP_VERSION", "v2.3.2")
    kraken2db_home = os.environ.get("KRAKEN2DB_HOME", "/ps/phoenix/kraken2db/")
    md5 = get_md5(wgs_input)

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

    # Create folders for output figures
    workflow_id = generate_workflow_id(md5, "MycoSNP", ds)
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
            "Process": "MycoSNP",
            "Category": 'WGS Workflow',
            "Created on": datetime.now(),
            "Status": "Processing"
        }
    )

    try:
        report_path = os.path.join(outdir, 'pipeline_info/execution_report.html')
        cmd = f"nextflow run CDCgov/mycosnp-nf -profile docker --workflow {workflow} --input {wgs_input} --outdir {outdir}"

        if (add_sra_file is not None and add_sra_file != ""):
            cmd += f" --add_sra_file {add_sra_file}"

        if (add_vcf_file is not None and add_vcf_file != ""):
            cmd += f" --add_vcf_file {add_vcf_file}"

        if (workflow == "PRE_MYCOSNP"):
            report_path = os.path.join(outdir, 'combined/pre-mycosnp_summary/pre-mycosnp-summary.csv')
            cmd += f" --assembler {assembler}"
            cmd += f" --shovill_depth {shovill_depth}"
            if (genome_size is not None and genome_size != ""):
                cmd += f" --genome_size {genome_size}"
            cmd += f" --min_contig_cov {min_contig_cov}"
            cmd += f" --min_contig_len  {min_contig_len}"
            if (gambit_db is not None and gambit_db != ""):
                cmd += f" --gambit_db {gambit_db}"
            if (gambit_h5_dir is not None and gambit_h5_dir != ""):
                cmd += f" --gambit_h5_dir {gambit_h5_dir}"
            if (subtype_db is not None and subtype_db != ""):
                cmd += f" --subtype_db {subtype_db}"
        elif (workflow == "NFCORE_MYCOSNP"):
            if (fasta is not None and fasta != ""):
                cmd += f" --fasta {fasta}"
            if (ref_dir is not None and ref_dir != ""):
                cmd += f" --ref_dir {ref_dir}"
            if (ref_masked_fasta is not None and ref_masked_fasta != ""):
                cmd += f" --ref_masked_fasta {ref_masked_fasta}"
            if (ref_fai is not None and ref_fai != ""):
                cmd += f" --ref_fai {ref_fai}"
            if (ref_dict is not None and ref_dict != ""):
                cmd += f" --ref_dict {ref_dict}"
            if (ref_bwa is not None and ref_bwa != ""):
                cmd += f" --ref_bwa {ref_bwa}"
            if (snpeffconfig is not None and snpeffconfig != ""):
                cmd += f" --snpeffconfig {snpeffconfig}"
            if (species is not None and species != ""):
                cmd += f" --species {species}"
            if (positions is not None and positions != ""):
                cmd += f" --positions {positions}"
            if (genes is not None and genes != ""):
                cmd += f" --genes {genes}"
            if (exclude is not None and exclude != ""):
                cmd += f" --exclude {exclude}"

            if (not save_reference):
                cmd += f" --save_reference false"
            if (not save_alignment):
                cmd += f" --save_alignment false"
            if (not rapidnj):
                cmd += f" --rapidnj false"
            if (not fasttree):
                cmd += f" --fasttree false"
            if (snpeff):
                cmd += f" --snpeff true"
            if (iqtree):
                cmd += f" --iqtree true"
            else:
                cmd += f" --iqtree false"
            if (raxmlng):
                cmd += f" --raxmlng true"
            else:
                cmd += f" --raxmlng false"
            if (save_debug):
                cmd += f" --save_debug true"
            else:
                cmd += f" --save_debug false"
            if (skip_samples is not None and skip_samples != ""):
                cmd += f" --skip_samples {skip_samples}"
            if (skip_samples_file is not None and skip_samples_file != ""):
                cmd += f" --skip_samples_file {skip_samples_file}"
            if (skip_combined_analysis):
                cmd += f" --skip_combined_analysis true"
            else:
                cmd += f" --skip_combined_analysis false"
            if (skip_phylogeny):
                cmd += f" --skip_phylogeny true"
            else:
                cmd += f" --skip_phylogeny false"

            cmd += f" --sample_ploidy {sample_ploidy}"
            cmd += f" --coverage  {coverage}"
            if (gvcfs_filter is not None and gvcfs_filter != ""):
                cmd += f" --gvcfs_filter {gvcfs_filter}"
            if (gatkgenotypes_filter is not None and gatkgenotypes_filter != ""):
                cmd += f" --gatkgenotypes_filter {gatkgenotypes_filter}"
            cmd += f" --max_amb_samples {max_amb_samples}"
            cmd += f" --max_perc_amb_samples  {max_perc_amb_samples}"
            cmd += f" --min_depth  {min_depth}"
            if (not mask):
                cmd += f" --mask false"

        # Execute the Nextflow command
        s = subprocess.call(cmd, shell=True)

        if os.path.exists(report_path):
            results = {
                "output": outdir,
                "workflow_id": workflow_id,
                "md5": md5,
                "wf_results": report_path,
                # "figures":fig_path, 
                # "process_ids": process_ids
            }
            
            upsert_jobs(
                {
                    "job_id": job_id, 
                    "output": outdir,
                    # "process_ids": process_ids,
                    "workflow_id": workflow_id,
                    "results": results,
                    "figures": fig_path, 
                    "Completed on": datetime.now(),
                    "Status": "Success"
                }
            )
        else:
            detail = f'MycoSNP workflow is failed: Report not generated.'
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
        detail = f"MycoSNP workflow is failed: {e}"
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