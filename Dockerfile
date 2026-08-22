# bdb-qchipms: Quantitative Locus-Specific Chromatin Proteomics Pipeline Container
FROM condaforge/mambaforge:latest

LABEL maintainer="Himanshu Bhandary <2032ushimanshu@gmail.com>"
LABEL description="Reproducible container for bdb-qchipms multi-omic pipeline execution"

# Set working directory
WORKDIR /opt/bdb-qchipms

# Copy environment configuration
COPY environment.yml .

# Install dependencies into base environment
RUN mamba env update -n base -f environment.yml && \
    mamba clean --all -f -y

# Copy repository code
COPY . .

# Set default entrypoint to snakemake
ENTRYPOINT ["snakemake"]
CMD ["--help"]
