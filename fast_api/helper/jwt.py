# app/auth/bearer.py
# Based on example provide by fastapi-jwt project team
# https://github.com/testdrivenio/fastapi-jwt

from typing import Optional

from fastapi import Request, HTTPException
from fastapi.security import HTTPBearer, HTTPAuthorizationCredentials
from fastapi_jwt_auth import AuthJWT


class JWTBearer(HTTPBearer):
    def __init__(self, auto_error: bool = True):
        super(JWTBearer, self).__init__(auto_error=auto_error)

    async def __call__(self, request: Request):
        credentials: HTTPAuthorizationCredentials = await super(
            JWTBearer, self
        ).__call__(request)
        if credentials:
            if not credentials.scheme == "Bearer":
                raise HTTPException(
                    status_code=403, detail="Invalid authentication scheme."
                )
            if not self.verify_jwt(request, credentials.credentials):
                raise HTTPException(
                    status_code=403, detail="Invalid token or expired token."
                )
            return credentials.credentials
        else:
            raise HTTPException(status_code=403, detail="Invalid authorization code.")

    # noinspection PyBroadException
    @staticmethod
    def verify_jwt(request: Request, token) -> bool:
        is_token_valid: bool = False
        auth_jwt = AuthJWT(request, token)
        try:
            payload = auth_jwt.get_raw_jwt()
        except Exception:
            payload = None
        if payload:
            is_token_valid = True
        return is_token_valid
