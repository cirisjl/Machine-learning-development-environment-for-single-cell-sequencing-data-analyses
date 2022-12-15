import fastapi
from fastapi.responses import JSONResponse
from cel_worker import start,test
from classes import Request,Example,WorkflowTaskResult,WorkflowTask
from celery.result import AsyncResult


router = fastapi.APIRouter()

@router.get("/example")
async def example(request: Example):
    task1 = test.delay(request.temp)

    return {"task_id":task1.task_id,"task_status":"asd"}

@router.post("/home")
async def trigger(req:Request):
    task = start.delay(req.input)

    return {"task_id":task.task_id,"status":"Processing"}

#
@router.get("/result/{task_id}",response_model=WorkflowTaskResult,status_code=202,responses={202:
                    {'model':WorkflowTask,'description':'Accepted: Not Ready'}})
async def workflow(task_id):
    task=AsyncResult(task_id)
    if not task.ready():
        return JSONResponse(status_code=202,content={"task_id":str(task_id),"task_status":"Processing"})
    result=task.get()
    return {'task_id':task_id,"task_status":"Success",'outcome':str(result)}
