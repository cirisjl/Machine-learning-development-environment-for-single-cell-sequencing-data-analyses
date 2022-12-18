from pydantic import BaseModel

class Request(BaseModel):
    dataset1 : str
    dataset2 : str

class WorkflowTask(BaseModel):
    task_id:str
    task_status:str

class WorkflowTaskResult(BaseModel):
    task_id:str
    task_status:str
    outcome:str