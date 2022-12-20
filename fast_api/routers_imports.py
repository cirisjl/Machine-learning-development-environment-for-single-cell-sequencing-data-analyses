# from routers.home import home
from routers.download import download
from rpy2_services import seurat_apis

routers = [
    download,
    seurat_apis
]
