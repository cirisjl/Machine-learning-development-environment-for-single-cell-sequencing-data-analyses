FROM rocker/r-ver:4.2.1



# Install OS packages
RUN apt-get update && apt-get install -yq --no-install-recommends htop vim unzip git wget nano rsync apt-utils parallel sshpass lsb-release libhdf5-dev hdf5-tools libgsl0-dev libxml2-dev libboost-all-dev libssl-dev libhdf5-dev curl libsodium-dev libudunits2-dev libgdal-dev libgeos-dev libproj-dev build-essential xorg-dev libreadline-dev libc6-dev zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libcairo2-dev libpango1.0-dev tcl-dev tk-dev openjdk-8-jdk gfortran libncurses5-dev libncursesw5-dev procps texlive libv8-dev libgit2-dev default-libmysqlclient-dev libpq-dev libsasl2-dev libsqlite3-dev libssh2-1-dev unixodbc-dev libpng-dev libfftw3-dev libgsl-dev libhiredis-dev gdebi-core ant clang cmake debsigs dpkg-sig expect fakeroot gnupg1 libacl1-dev libattr1-dev libcap-dev libclang-dev libunwind8-dev libegl1-mesa libfuse2 libgl1-mesa-dev libgtk-3-0 libpam-dev libuser1-dev libxslt1-dev lsof ninja-build patchelf pkg-config psmisc rrdtool software-properties-common uuid-dev  python3-dev python3-pip python3-tk python3-distutils && apt-get clean && rm -rf /var/lib/apt/lists
# WORKDIR /tmp

# Install UMAP
RUN LLVM_CONFIG=/usr/lib/llvm-10/bin/llvm-config pip3 install llvmlite
RUN pip3 install numpy==1.22.4
RUN pip3 install umap-learn

# Install FIt-SNE
RUN git clone --branch v1.2.1 https://github.com/KlugerLab/FIt-SNE.git
RUN g++ -std=c++11 -O3 FIt-SNE/src/sptree.cpp FIt-SNE/src/tsne.cpp FIt-SNE/src/nbodyfft.cpp  -o bin/fast_tsne -pthread -lfftw3 -lm


# Install python packages
RUN pip install --upgrade --no-cache awscli==1.16.14 colorama==0.3.9 PyYAML==3.13 protobuf==3.19.4 rich==11.1 cutadapt magic-impute scipy biopython h5py matplotlib matplotlib_venn tqdm python-igraph==0.9.11 scprep torch networkx pandas==1.3.5 seaborn munkres feather-format rpy2 anndata2ri bioinfokit leidenalg

# Install Jupyter
RUN pip install --upgrade --no-cache jupyter jupyter-server-proxy jupyter-rsession-proxy rpy2
RUN Rscript -e 'install.packages("IRkernel")' && \
    Rscript -e 'IRkernel::installspec()' 

# Install SAVER
RUN Rscript -e 'devtools::install_github("mohuangx/SAVER")'

# Install scGNN
RUN Rscript -e 'devtools::install_github("BMEngineeR/scGNNLTMG")'

# Install SIMLR
RUN Rscript -e 'devtools::install_github("BatzoglouLabSU/SIMLR", ref = "master")'

# Install scanpy&&dropkick&&scvi-tools
RUN pip install --upgrade --no-cache 'scanpy[leiden]' anndata muon dropkick scvi-tools

# Install snakemake
RUN pip3 install git+https://github.com/snakemake/snakemake

CMD ["bash"]
