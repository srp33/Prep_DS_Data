FROM bioconductor/bioconductor_docker:RELEASE_3_23

# Install system libraries commonly needed by R/Bioconductor packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R/Bioconductor packages (not included in the base image)
RUN R -e "BiocManager::install('tidyverse', ask = FALSE, update = FALSE)"

ADD install_BrainArray_packages.R /project/
ADD datasets.xlsx /project/
RUN Rscript /project/install_BrainArray_packages.R

RUN R -e "BiocManager::install('SCAN.UPC', ask = FALSE, update = FALSE)"
RUN R -e "BiocManager::install('GEOquery', ask = FALSE, update = FALSE)"
RUN R -e "BiocManager::install('biomaRt', ask = FALSE, update = FALSE)"