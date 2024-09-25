import base64
import gzip
import json
import pandas as pd
from io import BytesIO

def gzip_str(to_gzip: str) -> str:
    out = BytesIO()
    with gzip.GzipFile(fileobj=out, mode='w') as f: 
        f.write(to_gzip.encode())
    return base64.b64encode(out.getvalue()).decode()


def ungzip_str(to_ungzip: str) -> str:
    compressed = base64.b64decode(to_ungzip) 
    with gzip.GzipFile(fileobj=BytesIO(compressed)) as f: 
        return f.read().decode()
    

def gzip_dict(to_gzip: dict) -> str:
    jsonStr = json.dumps(to_gzip)
    return gzip_str(jsonStr)


def ungzip_dict(to_ungzip: str) -> dict:
    jsonStr = ungzip_str(to_ungzip)
    return json.loads(jsonStr)


def gzip_df(to_gzip: pd.DataFrame) -> str:
    # dropping null value columns to avoid errors 
    to_gzip.dropna(inplace = True) 
    
    # converting to dict 
    dfDict = to_gzip.to_dict() 
    return gzip_dict(dfDict)


def ungzip_df(to_ungzip: str) -> pd.DataFrame:
    dfDict = ungzip_dict(to_ungzip)
    return pd.DataFrame.from_dict(dfDict)