# Installation

## conda/mamba (recommended)

You can install carnation by creating a new conda environment containing
all dependencies as listed in the file `requirements.yaml`. This is *recommended*.

```
mamba env create -p env --file requirements.yaml
```

Next, activate the environment and install using `R CMD build` and `install.packages`:

```
conda activate ./env
R CMD build carnation/
Rscript -e "install.packages('carnation_1.0.tar.gz', repos=NULL)"
```

## devtools

You can also install carnation using the `devtools` package. Here, we use `devtools` to figure
out the dependencies directly, mimicking a direct R installation.

First, install `devtools` in an existing R installation, e.g. in RStudio.

```
install.packages('devtools')
```

If conda/mamba is installed, you could also do this with a conda environment
that contains only `devtools`, e.g.

```
conda create -p env r-devtools
conda activate ./env
```

Next, open R and run `setRepositories()` to set CRAN and Bioconductor repositories. I usually set the first 5 options,
CRAN, BioC software, BioC annotation, BioC experiment, CRAN (extras) to get dependencies
from both CRAN and Bioconductor.

Then, from inside the carnation source directory, run `devtools::install()` to install carnation. Note: This is
a time-consuming step, especially if working with a fresh R installation, as a large number of dependencies will be installed.

# Getting started

## First run

First, load the `carnation` library.

```
library(carnation)
```

`carnation` uses the `kaleido` python module to allow interactive plots to be
saved to PDF. To enable this functionality, you need to run:

```
install_carnation()
```

By default, this installs the python modules `plotly` and `kaleido` into a python virtual
environment named `r-carnation`.

Finally, run the main app function.

```
run_carnation()
```

This will open up carnation in a browser window. Note that, if you're running carnation
from a remote ssh server, you might need to add some additional options. For example,
if you're running carnation on a remote server that is forwarding port `12345` via
ssh tunneling:

```
run_carnation(options=list(port=12345, launch.browser=FALSE))
```

This will start carnation on port 12345 and you will be able to access it via
a web browser from the URL: `http://127.0.0.1:12345`.

## Data setup

Organize your data in a directory structure that carnation can use, e.g. in `/carnation/data`.
For example, let's say you have two projects *project1* and *project2*.

- *project1* has 3 analyses, *main*, *subset* and *main-nooutlier*
- *project2* has 1 analysis, *default*

Then your carnation data directory can look like this:

```
/carnation/data/
  ├─ project1
  │  ├─ main.rds
  │  ├─ subset.rds
  │  └─ main-nooutlier.rds
  │
  └─ project2
     └─ default.rds
```

The first time you run carnation, you will be asked to
choose a data area. If you have set up your data as above, type in the
location `/carnation/data` in the data area input.

If everything worked as expected, then the *Available projects* menu should be populated with two options *project1* and *project2*.

- If *project1* is selected, *Available assays* will show three options, *main*, *subset* and *main-nooutlier*.
- If *project2* is selected, there will be one available assay - *default*.

Subsequently, you can add data areas using the Settings tab in the app.
