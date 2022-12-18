import fastapi
from fastapi.responses import JSONResponse
from classes import Request,WorkflowTask,WorkflowTaskResult
from cel_worker import start
from celery.result import AsyncResult

router = fastapi.APIRouter()

@router.post("/batch_correction")
async def trigger(request:Request):
    task = start.delay(request.dataset1,request.dataset2)
    return{"Task_id":task.task_id,"Task_Status":"Processing"}

@router.get("/batch_correction/{task_id}",response_model=WorkflowTaskResult,status_code=202,responses={202:
                    {'model':WorkflowTask,'description':'Accepted: Not Ready'}})
async def workflow(task_id):
    task = AsyncResult(task_id)
    if not task.ready():
        return JSONResponse(status_code=202, content={"task_id": str(task_id), "task_status": "Processing"})
    result = task.get()
    return {'task_id': task_id, "task_status": "Success", 'outcome': str(result)}
