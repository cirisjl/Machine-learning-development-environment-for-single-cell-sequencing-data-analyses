FROM rocker/r-ver:4.2.1


# ENV R_LIBS=/usr/local/lib/R/site-library

# ENV R_LIBS_USER=/usr/local/lib/R/site-library
ENV RETICULATE_MINICONDA_ENABLED=FALSE

# Install R and required libraries
RUN apt-get update && \
    apt-get install -y python3 python3-pip && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Set Python path for reticulate
ENV RETICULATE_PYTHON /usr/bin/python3

RUN echo "options(defaultPackages = c(getOption('defaultPackages'), 'utils', 'grDevices', 'graphics', 'stats', 'methods'))" >> /usr/lib/R/etc/Rprofile.site

# Set up default R packages
# RUN echo "options(defaultPackages = c(getOption('defaultPackages'), 'utils', 'grDevices', 'graphics', 'stats', 'methods'))" >> /usr/lib/R/etc/Rprofile.site

# ## Install base R packages
# RUN Rscript -e 'install.packages(c("BiocManager"))' && \
#    Rscript -e 'BiocManager::install()'


RUN Rscript -e 'install.packages(c("dplyr", "R.utils", "cli", "Seurat", "reticulate"))'

# Install Seurat from GitHub
# RUN Rscript -e 'remotes::install_github("satijalab/seurat", ref = "release/4.0")'

# Install Bioconductor packages
# RUN Rscript -e 'BiocManager::install(c("org.Hs.eg.db"))'

# RUN Rscript -e 'BiocManager::install(c("SingleCellExperiment" , "EnsDb.Hsapiens.v86", "scater", "AnnotationDbi", "org.Hs.eg.db", "org.Mm.eg.db"))'

# Install github packages
# RUN Rscript -e 'devtools::install_github("stuart-lab/signac")'

RUN mkdir /app

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

# Set the FLASK_APP environment variable
ENV FLASK_APP=main.py

# Expose port 5003 for the Flask application to listen on
EXPOSE 5003

CMD [ "python3", "-m" , "flask", "run", "--host=0.0.0.0","--port=5003"]
