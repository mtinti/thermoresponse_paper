# One image for every step except the GPU embedding step.
#
# Build for linux/amd64: the published random forests were grown on Linux
# x86-64, and their splits depend on the platform's floating-point details.
#
#   docker build --platform linux/amd64 -t thermoresponse-paper .
#   docker run --rm --platform linux/amd64 -v "$PWD":/work thermoresponse-paper make all

FROM --platform=linux/amd64 mambaorg/micromamba:1.5.10

# R and Bioconductor for the differential-expression steps.
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.12 make \
        r-base=4.4.3 bioconductor-edger=4.4.0 bioconductor-limma=3.62.1 \
        bioconductor-cqn=1.52.0 \
    && micromamba clean --all --yes

# Python packages from PyPI wheels, the builds used for the published figure.
COPY requirements.txt /tmp/requirements.txt
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN pip install --no-cache-dir -r /tmp/requirements.txt

WORKDIR /work
