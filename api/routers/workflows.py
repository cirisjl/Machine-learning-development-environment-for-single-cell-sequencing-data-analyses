from celery import group
from fastapi import APIRouter
from starlette.responses import JSONResponse

# from api import tools
from celery_tasks.tasks import create_clustering_task
from schemas.schemas import Dataset
router = APIRouter(prefix='/api/workflows', tags=['workflows'], responses={404: {"description": "API Not found"}})


@router.post("/clustering")
async def create_clustering_task_async(ds: Dataset):
    """
    Create a task for clustering
    """
    ds_dict = ds.dict()  # Convert the Pydantic model to a dict
    task = create_clustering_task.apply_async(args=[ds_dict])
    return JSONResponse({"job_id": task.id, "status": "Clustering task submitted successfully"})
