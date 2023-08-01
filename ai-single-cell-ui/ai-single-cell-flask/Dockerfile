FROM python:3.9-slim-buster

# Install R and required libraries from the 'rocker/r-ver' image
RUN apt-get update && \
    apt-get install -y --no-install-recommends dirmngr gnupg && \
    echo "deb https://cloud.r-project.org/bin/linux/debian buster-cran40/" >> /etc/apt/sources.list && \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9 && \
    apt-get update && \
    apt-get install -y --no-install-recommends r-base r-base-dev libxml2-dev libssl-dev libcurl4-openssl-dev libcairo2-dev libxt-dev libbz2-dev liblzma-dev libffi-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Set up default R packages
RUN echo "options(defaultPackages = c(getOption('defaultPackages'), 'utils', 'grDevices', 'graphics', 'stats', 'methods', 'Seurat'))" >> /etc/R/Rprofile.site

## Install base R packages
RUN Rscript -e 'install.packages(c("BiocManager"))' && \
   Rscript -e 'BiocManager::install()'

RUN Rscript -e 'install.packages(c("dplyr", "R.utils", "reticulate", "Seurat","SeuratDisk","SeuratData", "patchwork", "cli"))'

# Install Bioconductor packages
# RUN Rscript -e 'BiocManager::install(c("org.Hs.eg.db"))'

RUN Rscript -e 'BiocManager::install(c("SingleCellExperiment" , "EnsDb.Hsapiens.v86", "scater", "AnnotationDbi", "org.Hs.eg.db", "org.Mm.eg.db"))'

# Set the default library location to /usr/local/lib/R/site-library
#RUN echo ".libPaths(c('/usr/local/lib/R/site-library', .libPaths()))" > /usr/local/lib/R/etc/Rprofile.site

# Install github packages
# RUN Rscript -e 'devtools::install_github("stuart-lab/signac")'

WORKDIR /usr/src/app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

# Set the FLASK_APP environment variable
ENV FLASK_APP=main.py

# Expose port 5003 for the Flask application to listen on
EXPOSE 5003

CMD [ "python3", "-m" , "flask", "run", "--host=0.0.0.0","--port=5003"]
