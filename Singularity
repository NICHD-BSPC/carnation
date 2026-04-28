Bootstrap: docker
From: condaforge/miniforge3:latest

%files
    requirements.yaml /app/carnation/

%environment
    export PATH=/env/carnation-env/bin:/opt/conda/bin/:$PATH
    export RETICULATE_PYTHON=/env/carnation-env/bin/python

%post -c /bin/bash
    /opt/conda/bin/conda env create -p /env/carnation-env --file /app/carnation/requirements.yaml
    /env/carnation-env/bin/Rscript -e "setRepositories(ind=1:5); remotes::install_github('NICHD-BSPC/carnation@r4.3', upgrade='never')"

%runscript
    #!/bin/bash

    exec /env/carnation-env/bin/R -e "library(carnation); run_carnation(options=list(launch.browser=FALSE, port=as.integer(Sys.getenv('PORT1'))))"
