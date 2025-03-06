FROM ubuntu:latest

RUN DEBIAN_FRONTEND=noninteractive apt-get update && apt-get install -y wget

WORKDIR /home/conda

# set shell to bash
SHELL ["/bin/bash", "-c"]

ENV SHELL /bin/bash

ENV CONDA_DIR /opt/conda

# Install miniforge
RUN wget --quiet https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O miniforge.sh && \
    /bin/bash miniforge.sh -b -p /opt/conda

# Put conda in path so we can use conda install
ENV PATH=$CONDA_DIR/bin:$PATH

# get carnation code
COPY . /app/carnation/

# build conda env
RUN conda env create -p /env/carnation-env --file /app/carnation/requirements-pinned.yaml

# set workdir to carnation
WORKDIR /app/carnation

# add conda env bin to path
ENV PATH=/env/carnation-env/bin:$PATH

# add path to python for reticulate
ENV RETICULATE_PYTHON=/env/carnation-env/bin/python

# add port to run carnation as command-line argument
ARG PORT=8080

# command to run carnation
CMD R -e "devtools::load_all(); run_carnation(options=list(host='127.0.0.1', port=as.integer(Sys.getenv('PORT'))))"

