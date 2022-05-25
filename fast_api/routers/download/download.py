from fastapi import APIRouter, Depends, HTTPException
from fastapi_jwt_auth import AuthJWT

from fastapi import Depends, HTTPException, status

from sqlalchemy.orm import Session

# import db.db_funcs as db
import helper.jwt as jwt
import helper.downloader_funcs as downloader
from pydantic import AnyUrl
from auth.credentials import PATH

router = APIRouter(
    prefix="/download",
    tags=["download"],
    # dependencies=[Depends(jwt.JWTBearer())],
    responses={404: {"description": "Not found"}},
)


@router.post(
    "/save",
)
async def post_download(
    url: AnyUrl,
):

    try:
        tmp = PATH + "tmp/"
        # path = "/Users/jsaied/GitHub/Machine-learning-development-environment-for-single-cell-sequencing-data-analyses/fast_api/tmp"
        download_path = await downloader.async_download(url=url, folder=tmp, log=True)
        # download_path = await downloader.download_url(url=url,folder=path, log=True)
        return {"message": "Downloaded successfully", "path": download_path}
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
