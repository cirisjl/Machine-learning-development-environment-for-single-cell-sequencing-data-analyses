import time
from fastapi import FastAPI, Request
from fastapi.middleware.cors import CORSMiddleware

# import db.db_funcs as db

from routers_imports import routers
from fastapi.responses import JSONResponse
from fastapi_jwt_auth.exceptions import AuthJWTException

from fastapi_jwt_auth import AuthJWT
from auth.credentials import ADMIN_USER, ADMIN_PASSWORD, PUB_KEY, PRIV_KEY, MODE
from pydantic import BaseModel

from fastapi.security import HTTPBasic
import helper.crypto as crypto
from fastapi import Depends, FastAPI, HTTPException, status
from fastapi.security import HTTPBasic, HTTPBasicCredentials

from fastapi.openapi.docs import get_swagger_ui_html
from fastapi.openapi.utils import get_openapi
from helper.openapi_help import custom_openapi
from fastapi.responses import RedirectResponse
import uvicorn
import jwt

from fastapi_jwt_auth import AuthJWT

#Db connection and related imports
from sqlalchemy.orm import Session

from db import crud, models, schemas, database
from db.database import SessionLocal, engine

security = HTTPBasic()


def login_validation(credentials: HTTPBasicCredentials = Depends(security)):

    if not (
        credentials.username == ADMIN_USER
        and crypto.verify_password(credentials.password, ADMIN_PASSWORD)
    ):
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Incorrect username or password",
            headers={"WWW-Authenticate": "Basic"},
        )
    return credentials.username

models.Base.metadata.create_all(bind=engine)  #Creating db session engine

app = FastAPI(
    title="FastAPI",
    description="Not Tested",
    version="1.0",
    docs_url=None,
    redoc_url=None,
    openapi_url=None,
)

# Dependency for db connection
def get_db():
    dbc = SessionLocal()
    try:
        yield dbc
    finally:
        dbc.close()


# @app.on_event("startup")
# async def startup():
#     db.async_engine = db.get_async_engine()


# @app.on_event("shutdown")
# async def shutdown_engine():
#     await db.async_engine.dispose()


app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


class Settings(BaseModel):
    authjwt_algorithm: str = "RS512"
    authjwt_public_key: str = PUB_KEY
    authjwt_private_key: str = PRIV_KEY


@AuthJWT.load_config
def get_config():
    return Settings()


for route in routers:
    app.include_router(route.router)


app.openapi = custom_openapi


@app.middleware("http")
async def add_process_time_header(request: Request, call_next):
    start_time = time.time()
    response = await call_next(request)
    process_time = time.time() - start_time
    response.headers["X-Process-Time"] = str(process_time)

    try:
        if request.headers.get("Authorization"):

            bearer_token = request.headers.get("Authorization")
            token_type = bearer_token.split(" ")[0]

            if not token_type == "Bearer":
                return response
            token = bearer_token.split(" ")[1]
            decoded = jwt.decode(token, options={"verify_signature": False})

            response.headers["X-User-Id"] = str(decoded["sub"])
    except Exception as e:
        print("Middleware Error: ", e)

    return response


@app.exception_handler(AuthJWTException)
def authjwt_exception_handler(request: Request, exc: AuthJWTException):
    return JSONResponse(status_code=exc.status_code, content={"detail": exc.message})


@app.get("/", include_in_schema=False)
async def get_redirect_to_swagger(username: str = Depends(login_validation)):
    return RedirectResponse(url="/docs")


@app.get("/docs", include_in_schema=False)
async def get_documentation(username: str = Depends(login_validation)):
    return get_swagger_ui_html(openapi_url="/openapi.json", title="docs")


@app.get("/openapi.json", include_in_schema=False)
async def openapi(username: str = Depends(login_validation)):
    return get_openapi(title="FastAPI", version="0.1.0", routes=app.routes)


#Demo APIs
@app.get("/getspatialdatasets")
async def getSpatialDatasets(db: Session = Depends(get_db)):
    all_datasets=crud.get_spatial_datasets(db)
    return all_datasets

@app.get("/getspatialdatasets_by_id/{id}")
async def getSpatialDatasets(id: str, db: Session = Depends(get_db)):
    dataset=crud.get_spatial_dataset_by_id(db,id)
    return dataset

@app.get("/getdatasetfactors")
async def getDatasetFactors(db: Session = Depends(get_db)):
    all_datasetfactors=crud.get_dataset_factors(db)
    return all_datasetfactors

@app.get("/getclusteringpublications")
async def getClusteringPublications(db: Session = Depends(get_db)):
    all_clustering_publications=crud.get_clustering_publications(db)
    return all_clustering_publications

@app.get("/getclusteringdatasetfactor")
async def getClusteringDatasetFactor(db: Session = Depends(get_db)):
    all_clustering_dataset_factor=crud.get_clustering_dataset_factor(db)
    return all_clustering_dataset_factor

if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=5006, debug=True)
