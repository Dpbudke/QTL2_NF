# QTL2_NF Pipeline - Comprehensive Container
# Single container for all r/qtl2 analysis modules
FROM rocker/r-ver:4.3.0

# Metadata
LABEL maintainer="Dpbudke <dpbudke@ucdavis.edu>"
LABEL description="Comprehensive container for QTL2_NF multiparental mouse QTL pipeline"
LABEL version="1.0.0"

# Install system dependencies for R packages and additional tools
RUN apt-get update && apt-get install -y \
    # Core R package dependencies
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
    # Additional system tools
    wget \
    curl \
    unzip \
    git \
    procps \
    # Cairo for high-quality graphics
    libcairo2-dev \
    libxt-dev \
    # For potential future genomics tools
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Install comprehensive R packages for QTL analysis
RUN R -e "install.packages(c( \
    # Core data manipulation and visualization \
    'dplyr', \
    'tidyr', \
    'readr', \
    'ggplot2', \
    'gplots', \
    'RColorBrewer', \
    'viridis', \
    # Statistical analysis \
    'qtl2', \
    'qtl2convert', \
    'qtl2plot', \
    'qtl2fst', \
    # Data quality and diagnostics \
    'corrplot', \
    'pheatmap', \
    'VennDiagram', \
    # Parallel processing \
    'parallel', \
    'foreach', \
    'doParallel', \
    # Utility packages \
    'devtools', \
    'here', \
    'optparse', \
    'yaml', \
    'jsonlite', \
    'data.table', \
    # Reporting \
    'knitr', \
    'rmarkdown' \
    ), repos='https://cran.rstudio.com/', dependencies=TRUE)"

# Install Bioconductor packages if needed for future modules
RUN R -e "if (!require('BiocManager', quietly = TRUE)) install.packages('BiocManager'); \
    BiocManager::install(c( \
    'GenomeInfoDb', \
    'rtracklayer', \
    'VariantAnnotation' \
    ), ask=FALSE, update=FALSE)"

# Create directory structure for pipeline scripts
RUN mkdir -p /opt/bin \
    && mkdir -p /opt/scripts \
    && mkdir -p /opt/data \
    && mkdir -p /data

# Copy all custom QTL pipeline scripts
# Module 1: Phenotype processing
COPY bin/robustZmat.R /opt/bin/robustZmat.R  
COPY bin/covCheck.R /opt/bin/covCheck.R

# Future modules (create placeholders or add as you develop)
# COPY bin/genoQC.R /opt/bin/genoQC.R
# COPY bin/crossValidation.R /opt/bin/crossValidation.R
# COPY bin/qtlPlots.R /opt/bin/qtlPlots.R

# Make all scripts executable and accessible
RUN chmod +x /opt/bin/*.R
ENV PATH="/opt/bin:${PATH}"

# Create a pipeline utilities script for common functions
RUN echo '#!/usr/bin/env Rscript' > /opt/bin/pipeline_utils.R && \
    echo '# Common utility functions for QTL2_NF pipeline' >> /opt/bin/pipeline_utils.R && \
    echo 'suppressPackageStartupMessages(library(qtl2))' >> /opt/bin/pipeline_utils.R && \
    echo 'suppressPackageStartupMessages(library(dplyr))' >> /opt/bin/pipeline_utils.R && \
    echo '' >> /opt/bin/pipeline_utils.R && \
    echo '# Function to validate sample IDs across files' >> /opt/bin/pipeline_utils.R && \
    echo 'validate_sample_ids <- function(pheno_ids, geno_ids = NULL) {' >> /opt/bin/pipeline_utils.R && \
    echo '  # Add validation logic here' >> /opt/bin/pipeline_utils.R && \
    echo '  return(TRUE)' >> /opt/bin/pipeline_utils.R && \
    echo '}' >> /opt/bin/pipeline_utils.R && \
    chmod +x /opt/bin/pipeline_utils.R

# Set R library path to ensure packages are found
ENV R_LIBS_USER="/usr/local/lib/R/site-library"

# Configure R for better performance in containers
RUN echo 'options(repos = c(CRAN = "https://cran.rstudio.com/"))' >> /usr/local/lib/R/etc/Rprofile.site && \
    echo 'options(download.file.method = "libcurl")' >> /usr/local/lib/R/etc/Rprofile.site && \
    echo 'options(Ncpus = parallel::detectCores())' >> /usr/local/lib/R/etc/Rprofile.site

# Set working directory
WORKDIR /data

# Health check to ensure R and key packages are working
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD R -e "library(qtl2); library(dplyr); cat('Container healthy\n')" || exit 1

# Default command - keep container running
CMD ["R", "--slave", "--no-restore", "--no-save"]

# Add build info as environment variables
ENV CONTAINER_BUILD_DATE="2024-12-09"
ENV CONTAINER_VERSION="1.0.0"
ENV R_VERSION="4.3.0"