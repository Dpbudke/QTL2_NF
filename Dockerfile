# QTL2_NF Pipeline - Comprehensive Container
# Single container for all r/qtl2 analysis modules
#
# Pinned so that rebuilding gives the R packages of image 1.0.0 (the lab's 2025 and 2026 runs; built April 6, 2026):
# the base image by digest, CRAN from the Posit Package Manager snapshot of April 7, 2026 (source packages, as in the
# original build; TIMBR's and Bioconductor's dependencies come from the base image's snapshot, as they did),
# Bioconductor 3.18 (the R 4.3 release) and TIMBR at commit 4d655aa. The last step compares every
# installed package with image_r_packages_1.0.0.tsv and fails on any difference. Ubuntu system libraries come from
# the 22.04 updates of the build day, so a rebuild is used only after a pilot gene set gives the same LOD scores as the
# image it replaces. Change software by editing this file and the list together, as a new image version.
FROM rocker/r-ver:4.3.0@sha256:48fb09f63e1cbcc1b0ce3974a8f206bff0804b6921bb36dfa08eafa264dad542

# Metadata
LABEL maintainer="Dpbudke <dpbudke@ucdavis.edu>"
LABEL description="Comprehensive container for QTL2_NF multiparental mouse QTL pipeline"
LABEL version="1.0.1"

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
    # GNU Scientific Library (required by TIMBR conditional posterior)
    libgsl-dev \
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
    'lineup2', \
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
    ), repos='https://packagemanager.posit.co/cran/2026-04-07', dependencies=TRUE)"

# Module 10: TIMBR allelic series analysis
RUN R -e "install.packages(c('ape', 'remotes'), repos='https://packagemanager.posit.co/cran/2026-04-07', dependencies=TRUE)"
# TIMBR's dependencies come from the base image's own package snapshot, as in image 1.0.0 (gmp, gsl, partitions)
RUN R -e "remotes::install_github('wesleycrouse/TIMBR@4d655aa0c0039a9928d73a190e4a0741b3bfee08', dependencies=TRUE, upgrade='never')"
RUN R -e "library(TIMBR); cat('TIMBR version:', as.character(packageVersion('TIMBR')), '\n')"

# Install Bioconductor packages if needed for future modules
RUN R -e "if (!require('BiocManager', quietly = TRUE)) install.packages('BiocManager'); \
    BiocManager::install(version = '3.18', ask = FALSE, update = FALSE); \
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

# Make scripts directory accessible (scripts are embedded in Nextflow process blocks)
RUN true
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
RUN echo 'options(repos = c(CRAN = "https://packagemanager.posit.co/cran/2026-04-07"))' >> /usr/local/lib/R/etc/Rprofile.site && \
    echo 'options(download.file.method = "libcurl")' >> /usr/local/lib/R/etc/Rprofile.site && \
    echo 'options(Ncpus = parallel::detectCores())' >> /usr/local/lib/R/etc/Rprofile.site

# Every package must have the version it had in image 1.0.0
COPY image_r_packages_1.0.0.tsv /opt/image_r_packages_1.0.0.tsv
RUN Rscript -e "want <- read.delim('/opt/image_r_packages_1.0.0.tsv', na.strings = 'NA'); \
    have <- installed.packages(); \
    got <- have[match(want\$Package, have[, 'Package']), 'Version']; \
    bad <- want\$Package[is.na(got) | got != want\$Version]; \
    if (length(bad)) stop('differs from image 1.0.0: ', paste(sprintf('%s %s (want %s)', bad, got[match(bad, want\$Package)], want\$Version[match(bad, want\$Package)]), collapse = ', ')); \
    cat(nrow(want), 'packages match image 1.0.0\n')"

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