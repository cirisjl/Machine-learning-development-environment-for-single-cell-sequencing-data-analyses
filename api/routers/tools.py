import os
from celery import group
from fastapi import APIRouter, HTTPException
from starlette.responses import JSONResponse

# from api import tools
from celery_tasks.tasks import create_qc_task, create_normalization_task, create_imputation_task, create_integration_task, create_evaluation_task, create_reduction_task, create_conversion_task, load_metadata_task, create_annotation_task
from schemas.schemas import Dataset, Datasets, PathRequest, UMAPRequest, UploadRequest
from datetime import datetime
from utils.mongodb import upsert_jobs

router = APIRouter(prefix='/api/tools', tags=['tools'], responses={404: {"description": "API Not found"}})


# @router.post("/ConvertToAnndata")
# async def ConvertToAnndata_task_async(request_data: PathRequest):
#     """
#     Convert Seurat/Single-Cell Experiment object to Anndata object and return the path of Anndata object or the list of assay names of Seurat object
#     """
#     path = request_data.path
#     adata_path,assay_names  = ConvertToAnndata_task(path)
#     if assay_names is None:
#         assay_names = []
#     if adata_path is None:
#         adata_path = "Not available"
#     # print("router")
#     # print("AssayNames")
#     # print(assay_names)
#     # print("adata_path")
#     # print(adata_path)

#     return JSONResponse({"assay_names": assay_names,"adata_path": adata_path , "message" : "OK"})

# @router.post("/qc")
# async def create_qc_task_async(ds: Dataset):
#     """
#     Create a task for quality control
#     """
#     task = create_qc_task.apply_async(args=[ds['dataset'], ds['input'], ds['userID'], ds['output'], ds['methods']], kwargs={'path_of_scrublet_calls':ds.path_of_scrublet_calls, 'show_error': ds.show_error})
#     return JSONResponse({"job_id": task.id})

def create_job(job_id, ds_dict: dict):
    upsert_jobs(
        {
            "job_id": job_id, 
            "description": ds_dict['description'],
            "datasetId": ds_dict['datasetId'],
            "method": ds_dict['method'],
            "datasetURL": ds_dict['input'],
            "process": ds_dict['process'],
            "class": 'category',
            "created_on": datetime.now(), 
            "status": "Queued"
        }
    )


@router.post("/qc")
async def create_qc_task_async(ds: Dataset):
    """
    Create a task for quality control
    """
    ds_dict = ds.dict()  # Convert the Pydantic model to a dict
    task = create_qc_task.apply_async(args=[ds_dict])
    create_job(task.id, ds_dict)

    return JSONResponse({"job_id": task.id, "status": "Quality Control task submitted successfully"})


@router.post("/normalize")
async def create_normalization_task_async(ds: Dataset):
    """
    Create a task for normalization
    """
    ds_dict = ds.dict()  # Convert the Pydantic model to a dict
    task = create_normalization_task.apply_async(args=[ds_dict])
    create_job(task.id, ds_dict)

    return JSONResponse({"job_id": task.id, "status": "Normalization task submitted successfully"})


@router.post("/impute")
async def create_imputation_task_async(ds: Dataset):
    """
    Create a task for imputation
    """
    ds_dict = ds.dict()  # Convert the Pydantic model to a dict
    task = create_imputation_task.apply_async(args=[ds_dict])
    create_job(task.id, ds_dict)

    return JSONResponse({"job_id": task.id, "status": "Imputation task submitted successfully"})


@router.post("/reduce")
async def create_reduction_task_async(ds: Dataset):
    """
    Create a task for imputation
    """
    ds_dict = ds.dict()  # Convert the Pydantic model to a dict
    task = create_reduction_task.apply_async(args=[ds_dict])
    create_job(task.id, ds_dict)

    return JSONResponse({"job_id": task.id, "status": "Dimension reduction task submitted successfully"})


@router.post("/convert")
async def create_conversion_task_async(ds: Dataset):
    """
    Create a task for imputation
    """
    ds_dict = ds.dict()  # Convert the Pydantic model to a dict
    task = create_conversion_task.apply_async(args=[ds_dict])
    create_job(task.id, ds_dict)

    return JSONResponse({"job_id": task.id, "status": "Conversion task submitted successfully"})


@router.post("/integrate")
async def create_integration_task_async(ids: Datasets):
    """
    Create a task for integration
    """
    ids_dict = ids.dict()  # Convert the Pydantic model to a dict
    task = create_integration_task.apply_async(args=[ids_dict])

    upsert_jobs(
        {
            "job_id": task.id, 
            "description": ids_dict['description'],
            "datasetIds": ids_dict['datasetIds'],
            "method": ids_dict['method'],
            "datasetURL": ids_dict['input'],
            "process": ids_dict['process'],
            "class": 'tools',
            "created_on": datetime.now(), 
            "status": "Queued"
        }
    )

    return JSONResponse({"job_id": task.id})


@router.post("/annotate")
async def create_annotation_task_async(ds: Dataset):
    """
    Create a task for imputation
    """
    ds_dict = ds.dict()  # Convert the Pydantic model to a dict
    task = create_annotation_task.apply_async(args=[ds_dict])
    create_job(task.id, ds_dict)

    return JSONResponse({"job_id": task.id, "status": "Annotation task submitted successfully"})


@router.post("/evaluate")
async def create_evaluation_task_async(ds: Dataset):
    """
    Create a task for evaluation
    """
    ds_dict = ds.dict() 
    task = create_evaluation_task.apply_async(args=[ds['dataset'], ds['input'], ds['userID'], ds['output'], ds['methods']], kwargs={'layer':ds['layer'], 'genes':ds.genes, 'ncores':ds.ncores, 'show_error': ds.show_error})
    create_job(task.id, ds_dict)

    return JSONResponse({"job_id": task.id})


@router.post('/metadata')
async def load_metadata_async(request_body: UploadRequest):
    file_dict = request_body.dict()
    task = load_metadata_task.apply_async(args=[file_dict])

    return JSONResponse({"job_id": task.id}) 