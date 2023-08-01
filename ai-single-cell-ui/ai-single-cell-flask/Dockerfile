FROM rocker/r-ver:4.2.1

# versions of software
ARG r_cran_version="cran40"
ARG rstudio_version="2023.03.0-386"

# Set global R options
RUN echo "options(repos = 'https://cloud.r-project.org')" > $(R --no-echo --no-save -e "cat(Sys.getenv('R_HOME'))")/etc/Rprofile.site
ENV RETICULATE_MINICONDA_ENABLED=FALSE

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

# WORKDIR /tmp

# Install RStudio
RUN RSTUDIO_PKG=rstudio-server-${rstudio_version}-amd64.deb && \
    cd /tmp && \
    wget -q https://download2.rstudio.org/server/bionic/amd64/${RSTUDIO_PKG} && \
    gdebi -n /tmp/${RSTUDIO_PKG} && \
    rm /tmp/${RSTUDIO_PKG}
ENV PATH="${PATH}:/usr/lib/rstudio-server/bin"
ENV LD_LIBRARY_PATH="/usr/lib/R/lib:/lib:/usr/lib/x86_64-linux-gnu:/usr/lib/jvm/java-8-openjdk-amd64/jre/lib/amd64/server"
RUN pip3 install radian


# Install base R packages
RUN Rscript -e 'install.packages(c("languageserver", "devtools", "Rmagic", "BiocManager"))' && \
    Rscript -e 'BiocManager::install()'

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
    "tinytex", "lsa", "uwot", "optparse", "DrImpute", "alluvial", "dplyr", "RColorBrewer"))'
# RUN Rscript -e 'tinytex::install_tinytex()'

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
ENV FLASK_APP=main.py``

# Expose port 5003 for the Flask application to listen on
EXPOSE 5003

CMD [ "python3", "-m" , "flask", "run", "--host=0.0.0.0","--port=5003"]
