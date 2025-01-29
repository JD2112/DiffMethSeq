FROM jd21/nano-dnam:latest
LABEL authors="Jyotirmoy Das" \
      description="Docker image containing all env reqs for edgeR and methylkit"

# build packages
ARG BUILD_PACKAGES="wget procps git pip"
ENV DEBIAN_FRONTEND noninteractive

RUN apt-get -q update && \
    apt-get -q install --yes ${BUILD_PACKAGES}

# install R packages required 
# copy and install the bioconductor packages into the image
COPY ./rpackages.R /rpackages.R
RUN R -e "install.packages('BiocManager', repos = 'http://cran.us.r-project.org')"
RUN Rscript /rpackages.R


ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8
