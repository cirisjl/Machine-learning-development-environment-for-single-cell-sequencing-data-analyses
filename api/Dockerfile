# FROM u15ec104/ai-single-cell:celery-backend-base
FROM stormgod/oscb-api

# Install python packages
RUN pip install --upgrade --no-cache \
    dash~=2.9.1 pymongo dash-bootstrap-components dash_bootstrap_templates python-multipart numpyro==0.12.1 rarfile loguru
    # louvain plotly psutil nvidia-ml-py pymde
    # gputil

# Downgrade SeuratObject
RUN R --no-echo --no-restore --no-save -e 'remotes::install_version("SeuratObject", "4.1.4", repos = c("https://satijalab.r-universe.dev", getOption("repos")))'
RUN R --no-echo --no-restore --no-save -e 'remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))'
# RUN R --no-echo --no-restore --no-save -e "install.packages('SeuratObject')"

# RUN mkdir /app
WORKDIR /app
COPY . .

# CMD ["bash"]
EXPOSE 5000
EXPOSE 5555