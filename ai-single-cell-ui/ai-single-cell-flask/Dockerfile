FROM rocker/r-ver:4.2.1

# versions of software
ARG r_cran_version="cran40"
ARG rstudio_version="2023.03.0-386"

# Set global R options
RUN echo "options(repos = 'https://cloud.r-project.org')" > $(R --no-echo --no-save -e "cat(Sys.getenv('R_HOME'))")/etc/Rprofile.site
ENV RETICULATE_MINICONDA_ENABLED=FALSE
# Set LD_LIBRARY_PATH
ENV LD_LIBRARY_PATH=/usr/local/lib/R/lib:$LD_LIBRARY_PATH


# Install OS packages
RUN apt-get update && apt-get install -yq --no-install-recommends \
    # Basic packages 
    htop vim unzip git wget nano rsync apt-utils sudo \
    parallel sshpass lsb-release libhdf5-dev hdf5-tools \
    # R dependecies
    libgsl0-dev libxml2-dev libboost-all-dev libssl-dev libhdf5-dev curl libsodium-dev \
    libudunits2-dev libgdal-dev libgeos-dev libproj-dev build-essential xorg-dev \
    libreadline-dev libc6-dev zlib1g-dev libbz2-dev liblzma-dev \
    libcurl4-openssl-dev libcairo2-dev libpango1.0-dev tcl-dev tk-dev openjdk-8-jdk \
    gfortran libncurses5-dev libncursesw5-dev procps texlive libv8-dev libgit2-dev \
    default-libmysqlclient-dev libpq-dev libsasl2-dev libsqlite3-dev libssh2-1-dev \
    unixodbc-dev libpng-dev libfftw3-dev libgsl-dev libhiredis-dev \
    # RStudio dependencies - https://github.com/rstudio/rstudio/tree/master/dependencies/linux
    gdebi-core ant clang cmake debsigs dpkg-sig expect fakeroot gnupg1 libacl1-dev \
    libattr1-dev libcap-dev libclang-6.0-dev libclang-dev libunwind8-dev \
    libegl1-mesa libfuse2 libgl1-mesa-dev libgtk-3-0 libpam-dev \
    libuser1-dev libxslt1-dev lsof ninja-build patchelf pkg-config \
    psmisc rrdtool software-properties-common uuid-dev \
    # libigraph0-dev libclang-8-dev libclang-6.0-dev llvm-10 \
    # Python packages
    python3-dev python3-pip python3-tk python3-distutils && \	
    # clean up
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

## Install base R packages
RUN Rscript -e 'install.packages(c("BiocManager", "devtools","remotes"))' && \
   Rscript -e 'BiocManager::install()'

RUN Rscript -e 'install.packages(c("dplyr", "R.utils", "reticulate", "Seurat", "patchwork", "cli", "anndata"))'

RUN Rscript -e 'BiocManager::install(c("SingleCellExperiment" , "EnsDb.Hsapiens.v86", "scater", "AnnotationDbi", "org.Hs.eg.db", "org.Mm.eg.db"))'

# Install github packages
# RUN Rscript -e 'devtools::install_github("stuart-lab/signac", "satijalab/seurat-data")'

# RUN Rscript -e 'remotes::install_github("mojaveazure/seurat-disk")'


WORKDIR /usr/src/app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

# Set the FLASK_APP environment variable
ENV FLASK_APP=main.py

# Expose port 5003 for the Flask application to listen on
EXPOSE 5003

CMD [ "python3", "-m" , "flask", "run", "--host=0.0.0.0","--port=5003"]
