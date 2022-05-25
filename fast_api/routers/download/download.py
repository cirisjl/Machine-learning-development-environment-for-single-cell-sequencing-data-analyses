from fastapi import APIRouter, Depends, HTTPException
from fastapi_jwt_auth import AuthJWT

from fastapi import Depends, HTTPException, status

from sqlalchemy.orm import Session

# import db.db_funcs as db
import helper.jwt as jwt
import helper.downloader_funcs as downloader
from pydantic import AnyUrl
from auth.credentials import PATH
from fastapi.responses import FileResponse
import os

router = APIRouter(
    prefix="/download",
    tags=["download"],
    # dependencies=[Depends(jwt.JWTBearer())],
    responses={404: {"description": "Not found"}},
)


@router.post("/", response_class=FileResponse)
async def post_download(
    url: AnyUrl,
):

    try:
        tmp = PATH + "tmp/"
        download_path = await downloader.async_download(url=url, folder=tmp, log=True)

        if os.path.exists(download_path):
            filename = url.rpartition("/")[2]
            return FileResponse(path=download_path, filename=filename)

    except Exception as e:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
