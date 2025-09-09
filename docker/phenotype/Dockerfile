FROM rocker/r-ver:4.3.0

# Install system dependencies for R packages
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R packages required for Module 1
RUN R -e "install.packages(c( \
    'dplyr', \
    'gplots', \
    'qtl2', \
    'readr', \
    'ggplot2' \
    ), repos='https://cran.rstudio.com/', dependencies=TRUE)"

# Create directories for custom scripts
RUN mkdir -p /opt/bin

# Copy custom QC scripts
COPY bin/robustZmat.R /opt/bin/robustZmat.R  
COPY bin/covCheck.R /opt/bin/covCheck.R

# Make scripts accessible
ENV PATH="/opt/bin:${PATH}"

# Set working directory
WORKDIR /data

# Default command
CMD ["R"]
