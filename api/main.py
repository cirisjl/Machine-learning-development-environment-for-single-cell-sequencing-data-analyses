import time
import uvicorn as uvicorn
from fastapi import FastAPI, WebSocket, Request, Query
from celery.result import AsyncResult
import asyncio
from fastapi.responses import HTMLResponse

from config.celery_utils import create_celery
from routers import tools, conversion
from fastapi.middleware.cors import CORSMiddleware
from fastapi.middleware.wsgi import WSGIMiddleware
from dash_app.dashboard import app as dashboard
from utils.logger import *
# from dash_app.dashboard import is_valid_query_param, get_dash_layout


def create_app() -> FastAPI:
    current_app = FastAPI(title="Asynchronous tasks processing with Celery and RabbitMQ",
                          description="AI-Ready for Single-Cell Data Analyses Event driven architecture"
                                      "FastAPI Application with Celery and RabbitMQ",
                          version="1.0.0", )

    current_app.celery_app = create_celery()
    current_app.include_router(tools.router)
    current_app.include_router(conversion.router)
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


@app.websocket("/log/{user_id}")
async def websocket_endpoint_log(websocket: WebSocket, user_id: str) -> None:
    """WebSocket endpoint for client connections

    Args:
        websocket (WebSocket): WebSocket request from client.
    """
    await websocket.accept()

    try:
        while True:
            await asyncio.sleep(1)
            logs = await log_reader(user_id, 30)
            await websocket.send_text(logs)
    except Exception as e:
        print(e)
    finally:
        await websocket.close()


@app.get("/status")
def get_status():
    return {"status": "ok"}


if __name__ == "__main__":
    uvicorn.run("main:app", host='0.0.0.0', port=5000, reload=True)
