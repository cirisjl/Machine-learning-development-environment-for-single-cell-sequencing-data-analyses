from celery import group
from fastapi import APIRouter, HTTPException
from starlette.responses import JSONResponse

# from api import tools
from celery_tasks.tasks import create_download_dataset_task
from schemas.schemas import DownloadDataset

router = APIRouter(prefix='/api/oscb-cli', tags=['oscb-cli'], responses={404: {"description": "API Not found"}})


@router.post("/downloadDataset")
async def download_dataset_task(ds: DownloadDataset):
    """
    Download a dataset after processing user specified tool
    """
    ds_dict = ds.dict()  # Convert the Pydantic model to a dict
    task = create_download_dataset_task.apply_async(args=[ds_dict])

    return JSONResponse({"job_id": task.id, "status": "Download Dataset Task Submitted Successfully!"})