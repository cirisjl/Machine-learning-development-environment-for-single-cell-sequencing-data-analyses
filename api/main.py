import time
import uvicorn as uvicorn
from fastapi import FastAPI, WebSocket
from celery.result import AsyncResult
import asyncio

from config.celery_utils import create_celery
from routers import tools
from fastapi.middleware.cors import CORSMiddleware
from fastapi.middleware.wsgi import WSGIMiddleware
from dash_app.dashboard import app as dashboard1
from app2 import app as dashboard2


def create_app() -> FastAPI:
    current_app = FastAPI(title="Asynchronous tasks processing with Celery and RabbitMQ",
                          description="AI-Ready for Single-Cell Data Analyses Event driven architecture"
                                      "FastAPI Application with Celery and RabbitMQ",
                          version="1.0.0", )

    current_app.celery_app = create_celery()
    # current_app.include_router(tools.router)
    return current_app


app = create_app()

# Mount the Dash app as a sub-application in the FastAPI server
app.mount("/dashboard2", WSGIMiddleware(dashboard2.server))

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
    print('inside middleware!')
    start_time = time.time()
    response = await call_next(request)
    process_time = time.time() - start_time
    response.headers["X-Process-Time"] = str(f'{process_time:0.4f} sec')
    return response


@app.websocket("/taskStatus/{taskIdsCommaSeparated}")
async def websocket_endpoint(websocket: WebSocket, taskIdsCommaSeparated: str):
    await websocket.accept()
    while True:
        taskIds = taskIdsCommaSeparated.split(',')
        results = {}
        for taskId in taskIds:
            result = AsyncResult(taskId)
            if result.ready():
                if result.successful():
                    results[taskId] = 'Success'
                else:
                    results[taskId] = 'Failed'
            else:
                results[taskId] = 'Processing'
        await websocket.send_json(results)
        await asyncio.sleep(3)


@app.get("/status")
def get_status():
    return {"status": "ok"}


if __name__ == "__main__":
    uvicorn.run("main:app", host='0.0.0.0', port=5000, reload=True)
