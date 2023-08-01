FROM python:3.9-slim-buster

# Install R and required libraries from the 'rocker/r-ver' image
RUN apt-get update && \
    apt-get install -y --no-install-recommends dirmngr gnupg && \
    echo "deb https://cloud.r-project.org/bin/linux/debian buster-cran40/" >> /etc/apt/sources.list && \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9 B8F25A8A73EACF41 && \
    apt-get update && \
    apt-get install -y --no-install-recommends r-base r-base-dev libssh2-1-dev libgit2-dev libxml2-dev libssl-dev libcurl4-openssl-dev libcairo2-dev libxt-dev libbz2-dev liblzma-dev libffi-dev && \
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
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Set up default R packages
RUN echo "options(defaultPackages = c(getOption('defaultPackages'), 'utils', 'grDevices', 'graphics', 'stats', 'methods', 'Seurat'))" >> /etc/R/Rprofile.site

## Install base R packages
RUN Rscript -e 'install.packages(c("BiocManager", "devtools"))' && \
   Rscript -e 'BiocManager::install()'

# RUN Rscript -e 'install.packages(c("reticulate", "Seurat", "anndata", "patchwork"))'

# Install other CRAN
RUN Rscript -e 'install.packages(c( \
    "remotes", "RCurl", "Seurat", "rJava", "umap", "bookdown", "cluster", "KernSmooth", \
    "rmarkdown", "ROCR", "googleVis", "ggbeeswarm", "SLICER", "ggfortify", "VGAM", "feather", \
    "tidyverse", "pheatmap", "plyr", "dplyr", "readr", "reshape", "Rcpp", "Matrix", \
    "reshape2", "reticulate", "viridis", "ggseqlogo", "ggplot2", "ggthemes", "cowplot", \
    "ggforce", "ggridges", "ggrepel", "gplots", "igraph", "car", "Rfast2", "RcppAnnoy", \
    "ggpubr", "httpuv", "xtable", "sourcetools", "modeltools", "R.oo", "pracma", "knitr", \
    "R.methodsS3", "shiny", "later", "checkmate", "bibtex", "lsei", "ape", "RSpectra", \
    "bit", "segmented", "mclust", "flexmix", "prabclus", "diptest", "mvtnorm", \
    "robustbase", "kernlab", "trimcluster", "proxy", "R.utils", "htmlwidgets", \
    "hexbin", "crosstalk", "promises", "acepack", "zoo", "npsurv", "iterators", \
    "snow", "bit64", "permute", "mixtools", "lars", "ica", "fpc", "enrichR", "SDMTools", \
    "pbapply", "irlba", "dtw", "plotly", "metap", "lmtest", "fitdistrplus", "Hmisc", "png", \
    "foreach", "vegan", "tidyr", "withr", "magrittr", "rmpi", "knitr", "future", "kableExtra", \
    "statmod", "mvoutlier", "penalized", "mgcv", "corrplot", "scales", "rliger", "pandoc", \
    "rlang", "R6", "cli", "crayon", "SeuratObject", "stringi", "anndata", "loomR", \
    "tinytex", "lsa", "uwot", "optparse", "DrImpute", "alluvial", "dplyr", "RColorBrewer", "patchwork"))'

# Install Bioconductor packages
# RUN Rscript -e 'BiocManager::install(c("org.Hs.eg.db"))'

# RUN Rscript -e 'BiocManager::install(c("SingleCellExperiment" , "EnsDb.Hsapiens.v86", "scater", "AnnotationDbi", "org.Hs.eg.db", "org.Mm.eg.db"))'

# Set the default library location to /usr/local/lib/R/site-library
#RUN echo ".libPaths(c('/usr/local/lib/R/site-library', .libPaths()))" > /usr/local/lib/R/etc/Rprofile.site

# Install github packages
# RUN Rscript -e 'devtools::install_github("stuart-lab/signac", "mojaveazure/seurat-disk", "satijalab/seurat-data")'

# Install Bioconductor packages
RUN Rscript -e 'BiocManager::install(c( \
    "multtest", "S4Vectors","graph", "RBGL", "gtools", "xtable", "pcaMethods", "limma", "SingleCellExperiment", "Dict", \
    "Rhdf5lib", "GenomicRanges", "scater", "IRanges", "scran", "RUVSeq", "sva", "SC3", "TSCAN", "monocle", "destiny", \
    "rtracklayer", "BiocGenerics", "DESeq2", "edgeR", "sctransform", "MAST", "scmap", "biomaRt", "MultiAssayExperiment", "SummarizedExperiment", \
    "Biobase","monocle", "beachmat", "DropletUtils", "EnsDb.Hsapiens.v86", "EnsDb.Mmusculus.v79", "batchelor", "SingleR", "celldex", "singleCellTK", \
    "convert", "glmGamPoi", "AnnotationDbi", "org.Hs.eg.db", "org.Mm.eg.db"))'

# Install github packages
RUN Rscript -e 'devtools::install_github(c( \
    "immunogenomics/harmony", "tallulandrews/M3Drop", "hemberg-lab/scRNA.seq.funcs", \
    "Vivianstats/scImpute", "theislab/kBET", "kieranrcampbell/ouija", "hemberg-lab/scfind", \
    "cole-trapnell-lab/monocle3", "mojaveazure/seurat-disk", "satijalab/seurat-wrappers", \
    "satijalab/seurat-data", "stuart-lab/signac"))'

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
