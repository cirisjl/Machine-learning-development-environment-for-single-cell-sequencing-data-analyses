import time
import uvicorn as uvicorn
from fastapi import FastAPI, WebSocket, Request, Query
from celery.result import AsyncResult
import asyncio
from fastapi.responses import HTMLResponse

from config.celery_utils import create_celery
from routers import tools, benchmarks, worlkflows
from config.celery_utils import get_task_info
from fastapi.middleware.cors import CORSMiddleware
from fastapi.middleware.wsgi import WSGIMiddleware
from dash_app.dashboard import app as dashboard
from utils.redislogger import *
# from dash_app.dashboard import is_valid_query_param, get_dash_layout


def create_app() -> FastAPI:
    current_app = FastAPI(title="Asynchronous tasks processing with Celery and RabbitMQ",
                          description="AI-Ready for Single-Cell Data Analyses Event driven architecture"
                                      "FastAPI Application with Celery and RabbitMQ",
                          version="1.0.0", )

    current_app.celery_app = create_celery()
    current_app.include_router(tools.router)
    current_app.include_router(benchmarks.router)
    return current_app


app = create_app()

# Mount the Dash app as a sub-application in the FastAPI server
app.mount("/dashboard", WSGIMiddleware(dashboard.server))

celery = app.celery_app

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.middleware("http")
async def add_process_time_header(request, call_next):
    start_time = time.time()
    response = await call_next(request)
    process_time = time.time() - start_time
    response.headers["X-Process-Time"] = str(f'{process_time:0.4f} sec')
    return response


@app.websocket("/{request_type}/{job_id}")
async def websocket_endpoint(websocket: WebSocket, request_type:str, job_id: str):
    await websocket.accept()
    try:
        while True:
            if request_type == 'taskCurrentStatus':
                result = get_task_info(job_id)
                await websocket.send_json(result)
            elif request_type == 'log':
                logs = await log_reader(job_id, 0, 30)
                await websocket.send_text(logs)
            await asyncio.sleep(3)
    except Exception as e:
        print(e)
    finally:
        await websocket.close()

        
# @app.websocket("/{request_type}/{job_idsCommaSeparated}")
# async def websocket_endpoint(websocket: WebSocket, request_type:str, job_idsCommaSeparated: str):
#     await websocket.accept()
    
#     try:
#         while True:
#             if request_type == 'taskStatus':
#                 job_ids = job_idsCommaSeparated.split(',')
#                 results = {}
#                 for job_id in job_ids:
#                     result = AsyncResult(job_id)
#                     if result.ready():
#                         if result.successful():
#                             results[job_id] = 'Success'
#                         else:
#                             results[job_id] = 'Failed'
#                     else:
#                         results[job_id] = 'Processing'
#                 await websocket.send_json(results)

#             elif request_type == 'log':
#                 logs = await log_reader(job_idsCommaSeparated, 0, 30)
#                 await websocket.send_text(logs)
#             await asyncio.sleep(3)
#     except Exception as e:
#         print(e)
#     finally:
#         await websocket.close()


@app.get("/status")
def get_status():
    return {"status": "ok"}


@app.get("/api/task/{job_id}")
async def get_task_status(job_id: str) -> dict:
    """
    Return the status of the submitted Task
    """
    return get_task_info(job_id)


if __name__ == "__main__":
    uvicorn.run("main:app", host='0.0.0.0', port=5005, reload=True)
