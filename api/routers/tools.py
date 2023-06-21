from celery import group
from fastapi import APIRouter
from starlette.responses import JSONResponse

# from api import tools
from celery_tasks.tasks import create_qc_task, create_normalization_task, create_imputation_task
from config.celery_utils import get_task_info
from schemas.schemas import Dataset
router = APIRouter(prefix='/tools', tags=['tool'], responses={404: {"description": "Not found"}})


@router.post("/qc")
async def create_qc_task_async(ds: Dataset):
    """
    Create a task for quality control
    """
    task = create_qc_task.apply_async(args=[ds.dataset, ds.input, ds.userID, ds.output, ds.methods], kwargs={'path_of_scrublet_calls':ds.path_of_scrublet_calls, 'show_error': ds.show_error})
    return JSONResponse({"task_id": task.id})


@router.post("/normalize")
async def create_normalization_task_async(ds: Dataset):
    """
    Create a task for normalization
    """
    task = create_normalization_task.apply_async(args=[ds.dataset, ds.input, ds.userID, ds.output, ds.methods], kwargs={'default_assay':ds.default_assay, 'output_format':ds.output_format, 'species':ds.species, 'idtype':ds.idtype, 'show_umap': ds.show_umap, 'show_error': ds.show_error})
    return JSONResponse({"task_id": task.id})


@router.post("/impute")
async def create_imputation_task_async(ds: Dataset):
    """
    Create a task for imputation
    """
    task = create_imputation_task.apply_async(args=[ds.dataset, ds.input, ds.userID, ds.output, ds.methods], kwargs={'layer':ds.layer, 'genes':ds.genes, 'ncores':ds.ncores, 'show_error': ds.show_error})
    return JSONResponse({"task_id": task.id})


@router.get("/task/{task_id}")
async def get_task_status(task_id: str) -> dict:
    """
    Return the status of the submitted Task
    """
    return get_task_info(task_id)

