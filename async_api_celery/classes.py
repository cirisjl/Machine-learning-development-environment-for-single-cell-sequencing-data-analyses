from pydantic import BaseModel

class Request(BaseModel):
    input: str

class Example(BaseModel):
    temp : str

class WorkflowTask(BaseModel):
    task_id:str
    task_status:str

class WorkflowTaskResult(BaseModel):
    task_id:str
    task_status:str
    outcome:str