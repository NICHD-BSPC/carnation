FROM condaforge/miniforge3

# set shell to bash
SHELL ["/bin/bash", "-c"]

ENV SHELL=/bin/bash

ENV CONDA_DIR=/opt/conda

# Put conda in path so we can use conda install
ENV PATH=$CONDA_DIR/bin:$PATH

# get carnation code
COPY . /app/carnation/

# build conda env
RUN conda env create -p /env/carnation-env --file /app/carnation/requirements-pinned.yaml

# set workdir to carnation
WORKDIR /app

# add conda env bin to path
ENV PATH=/env/carnation-env/bin:$PATH

# set env variable for docker
ENV IN_DOCKER=true

# add path to python for reticulate
ENV RETICULATE_PYTHON=/env/carnation-env/bin/python

# add port to run carnation as command-line argument
ENV PORT=3838

EXPOSE 3838

# install carnation into env
RUN Rscript -e "setRepositories(ind=1:5); remotes::install_github('NICHD-BSPC/carnation', upgrade='never')"

# run carnation
CMD R -e "library(carnation); run_carnation(options=list(host='0.0.0.0', port=as.integer(Sys.getenv('PORT'))))"

