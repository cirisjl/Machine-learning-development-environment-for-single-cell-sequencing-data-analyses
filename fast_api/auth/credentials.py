import yaml
import os
from dotenv import load_dotenv

load_dotenv()

if os.getenv("place") == "local":
    path = "D:/Master's/ML_backend_AISingleCell/Machine-learning-development-environment-for-single-cell-sequencing-data-analyses/fast_api/"
elif os.getenv("place") == "docker":
    path = "/code/"
else:
    path = "D:/Master's/ML_backend_AISingleCell/Machine-learning-development-environment-for-single-cell-sequencing-data-analyses/fast_api/"


with open(path + "config.yaml") as file:
    yaml_file = yaml.load(file, Loader=yaml.FullLoader)
    db = yaml_file["database"]
    admin = yaml_file["admin"]
    mode = yaml_file["mode"]


db_credentials = db

ADMIN_USER = admin["user"]
ADMIN_PASSWORD = admin["hashed_password"]


PRIV_KEY = open(os.path.join(path, "auth", "RS512.key"), "rb").read().decode("utf-8")
PUB_KEY = open(os.path.join(path, "auth", "RS512.key.pub"), "rb").read().decode("utf-8")


BEARER_TOKEN = os.getenv("token", "")
SERVER_TOKEN = os.getenv("server_token", "")

MODE = mode
PATH = path
