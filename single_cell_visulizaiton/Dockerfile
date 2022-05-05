# Install Seurat 4.1.0
FROM satijalab/seurat:4.1.0

# Install Ubuntu packages
RUN apt-get update && apt-get install -y \
    gdebi-core \
    sudo \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libxt-dev \
    libssl-dev \
    libcairo2-dev \
	libsasl2-dev \
	libxml2-dev \
	libgit2-dev \
	default-libmysqlclient-dev \
	libpq-dev \
	libsqlite3-dev \
	libssh2-1-dev \
	unixodbc-dev \
	python \
	python3-pip
	
RUN apt-get -y --no-install-recommends install \
	libbz2-dev \
	liblzma-dev \
	libsodium-dev \
	libhiredis-dev

# Copy configuration files into the Docker image
COPY shiny-server.conf  /etc/init/shiny-server.conf

# Download and install ShinyServer (latest version)
RUN wget --no-verbose https://s3.amazonaws.com/rstudio-shiny-server-os-build/ubuntu-12.04/x86_64/VERSION -O "version.txt" && \
    VERSION=$(cat version.txt)  && \
    wget --no-verbose "https://s3.amazonaws.com/rstudio-shiny-server-os-build/ubuntu-12.04/x86_64/shiny-server-$VERSION-amd64.deb" -O ss-latest.deb && \
    gdebi -n ss-latest.deb && \
    rm -f version.txt ss-latest.deb

COPY /app /srv/shiny-server/

# Install R packages that are required
# TODO: add further package if you need!
RUN R -e "install.packages(c('seriation', 'dendextend', 'shinyjs', 'RColorBrewer', 'tidyr', 'tibble', 'dplyr', 'logging', 'shinyBS', 'periscope', 'shiny', 'shinydashboard', 'gdata', 'heatmaply', 'mongolite', 'plotly', 'Matrix', 'canvasXpress', 'future', 'pbapply'), repos=c('http://cran.rstudio.com/', 'https://github.com/jeroen/mongolite/'))"

# Make the ShinyApp available at port 3838
EXPOSE 3838

# Copy further configuration files into the Docker image
COPY shiny-server.sh /usr/bin/shiny-server.sh

RUN ["chmod", "+x", "/usr/bin/shiny-server.sh"]
RUN ["chmod", "777", "/srv/shiny-server/counter.Rdata"]
CMD ["/usr/bin/shiny-server.sh"]

ENV WAIT_VERSION 2.7.2
ADD https://github.com/ufoscout/docker-compose-wait/releases/download/$WAIT_VERSION/wait /wait
RUN chmod +x /wait
