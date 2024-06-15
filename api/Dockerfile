# FROM u15ec104/oscb:celery-backend-base
# FROM stormgod/oscb-api
# FROM stormgod/oscb-seurat:5.02
FROM stormgod/oscb-seurat:4.4

# RUN apt-get update && apt-get install -yq --no-install-recommends libglpk-dev libmagick++-dev

# Install python packages
# RUN pip install --upgrade --no-cache \
#     dash~=2.9.1 pymongo dash-bootstrap-components dash_bootstrap_templates python-multipart numpyro==0.12.1 rarfile loguru wsproto attrdict
    # louvain plotly psutil nvidia-ml-py pymde
    # gputil

# Downgrade SeuratObject
# RUN R --no-echo --no-restore --no-save -e 'remotes::install_version("SeuratObject", "4.1.4", repos = c("https://satijalab.r-universe.dev", getOption("repos")))'
# RUN R --no-echo --no-restore --no-save -e 'remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))'
# RUN R --no-echo --no-restore --no-save -e "install.packages('SeuratObject')"
# RUN R --no-echo --no-restore --no-save -e "install.packages('igraph')"
# RUN R --no-echo --no-restore --no-save -e "install.packages("irlba", type="source", force=TRUE)"
# RUN R --no-echo --no-restore --no-save -e "install.packages('magick')"
# RUN Rscript -e 'BiocManager::install("singleCellTK")'

# Install DoubletFinder
# RUN R --no-echo --no-restore --no-save -e "remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')"

# RUN pip install --upgrade --no-cache json-numpy

# RUN mkdir /app
# RUN pip install --upgrade --no-cache celery_once
# RUN R --no-echo --no-restore --no-save -e "install.packages('redux')"
# RUN R --no-echo --no-restore --no-save -e 'remotes::install_version("matrixStats", version="1.1.0")'
RUN pip install --upgrade --no-cache boltons

RUN R --no-echo --no-restore --no-save -e "install.packages('redux')"
WORKDIR /app
COPY . .

# CMD ["bash"]
EXPOSE 5000
EXPOSE 5555