from celery import group
from fastapi import APIRouter
from starlette.responses import JSONResponse

# from api import tools
from celery_tasks.tasks import create_qc_task, create_normalization_task, create_imputation_task, create_integration_task, create_evaluation_task, create_reduction_task
from config.celery_utils import get_task_info
from schemas.schemas import Dataset, IntegrationDataset, PathRequest
router = APIRouter(prefix='/tools', tags=['tool'], responses={404: {"description": "API Not found"}})


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
#     task = create_qc_task.apply_async(args=[ds.dataset, ds.input, ds.userID, ds.output, ds.methods], kwargs={'path_of_scrublet_calls':ds.path_of_scrublet_calls, 'show_error': ds.show_error})
#     return JSONResponse({"task_id": task.id})


@router.post("/qc")
async def create_qc_task_async(ds: Dataset):
    """
    Create a task for quality control
    """
    print("inside the API")
    ds_dict = ds.dict()  # Convert the Pydantic model to a dict
    task = create_qc_task.apply_async(args=[ds_dict])
    return JSONResponse({"task_id": task.id, "status": "Quality Control task submitted successfully"})


@router.post("/normalize")
async def create_normalization_task_async(ds: Dataset):
    """
    Create a task for normalization
    """
    ds_dict = ds.dict()  # Convert the Pydantic model to a dict
    task = create_normalization_task.apply_async(args=[ds_dict])
    return JSONResponse({"task_id": task.id, "status": "Normalization task submitted successfully"})


@router.post("/impute")
async def create_imputation_task_async(ds: Dataset):
    """
    Create a task for imputation
    """
    ds_dict = ds.dict()  # Convert the Pydantic model to a dict
    task = create_imputation_task.apply_async(args=[ds_dict])
    return JSONResponse({"task_id": task.id, "status": "Imputation task submitted successfully"})


@router.post("/impute")
async def create_reduction_task_async(ds: Dataset):
    """
    Create a task for imputation
    """
    ds_dict = ds.dict()  # Convert the Pydantic model to a dict
    task = create_imputation_task.apply_async(args=[ds_dict])
    return JSONResponse({"task_id": task.id, "status": "Dimension reduction task submitted successfully"})


@router.post("/integrate")
async def create_integration_task_async(ds: IntegrationDataset):
    """
    Create a task for integration
    """
    task = create_integration_task.apply_async(args=[ds.dataset, ds.input, ds.userID, ds.output, ds.methods, ds.species], kwargs={'default_assay':ds.default_assay, 'output_format':ds.output_format, 'genes':ds.genes, 'reference':ds.reference, 'show_error': ds.show_error})
    return JSONResponse({"task_id": task.id})


@router.post("/evaluate")
async def create_evaluation_task_async(ds: Dataset):
    """
    Create a task for evaluation
    """
    task = create_evaluation_task.apply_async(args=[ds.dataset, ds.input, ds.userID, ds.output, ds.methods], kwargs={'layer':ds.layer, 'genes':ds.genes, 'ncores':ds.ncores, 'show_error': ds.show_error})
    return JSONResponse({"task_id": task.id})


@router.get("/task/{task_id}")
async def get_task_status(task_id: str) -> dict:
    """
    Return the status of the submitted Task
    """
    return get_task_info(task_id)

