Bootstrap: docker
From: ubuntu:latest

%files
    ./R/* /app/carnation/R/
    ./inst /app/carnation
    ./man /app/carnation
    NAMESPACE /app/carnation
    DESCRIPTION /app/carnation
    .Rprofile /app/carnation
    env.yml /app/carnation

%environment
    export PATH=/env/carnation-env/bin:/opt/conda/bin/:$PATH

%post -c /bin/bash
    apt-get update && apt-get install -y wget
    wget --quiet https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O miniforge.sh && \
        /bin/bash miniforge.sh -b -p /opt/conda

    /opt/conda/bin/conda env create -p /env/carnation-env --file /app/carnation/env.yml

%runscript
    #!/bin/bash

    exec /env/carnation-env/bin/R -e "devtools::load_all('/app/carnation/'); run_carnation(options=list(launch.browser=FALSE, port=as.integer(Sys.getenv('PORT1'))))"
