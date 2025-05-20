# Carnation

Interactive shiny app for deeply exploring bulk RNA-Seq data, incorporating functional enrichment
and pattern analysis.

- Assess differential expression analysis results using interactive PCA, UpSet plots, heatmaps
  and a highly customizable gene plot.
- Dive into functional enrichment analysis with the help of fuzzy search across genes or pathway descriptions.
- Easily track genes of interest across the app using the novel "Gene scratchpad".
- Data access interface with an inbuilt authentication layer can be used in server mode to share results with
  others or to organize RNA-Seq projects on your local machine.

# Installation

## conda (recommended)

You can install carnation by creating a new conda environment containing
all dependencies as listed in the file `requirements-pinned.yaml`. This is *recommended*.

**WARNING:** Make sure you build the conda environment *outside* the carnation directory
or the build step will fail.

```
cd .. && conda env create -p env --file carnation/requirements-pinned.yaml
```

Next, activate the environment and start R:

```
conda activate ./env
R
```

Finally, install using `devtools::install_github`:

```
devtools::install_github('NICHD-BSPC/carnation')
```

## remotes

You can also install carnation using the `remotes` package. Here, we use `remotes` to figure
out the dependencies directly, mimicking a direct R installation.

First, install `remotes` in an existing R installation, e.g. in RStudio.

```
install.packages('remotes')
```

If conda is installed, you could also do this with a conda environment
that contains only `remotes`, e.g.

```
conda create -p env r-remotes
conda activate ./env
```

Next, open R, set repositories to get both CRAN and bioconductor packages and run
`remotes::install_github`.

```
setRepositories(ind=c(1,2,3,4,5))
remotes::install_github('NICHD-BSPC/carnation')
```

Note: This may be a time-consuming step, especially if working with a fresh R installation, as a large number of dependencies will be installed.

# Getting started

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

This will start carnation on port `12345`.

Now point your web browser to the URL: `http://127.0.0.1:12345`.
Since, this is your first run, carnation will ask you where your data is via
a modal dialog.


# Server mode

Carnation supports running in server mode using `shinymanager` to add an
authentication layer. For this, you need to first create a local sqlite database
with user details (for more details see [shinymanager docs](https://datastorm-open.github.io/shinymanager/)).

```
# create data frame with user details
credentials <- data.frame(
  user = c('shinymanager'),
  password = c('12345'),
  # password will automatically be hashed
  admin = c(TRUE),
  stringsAsFactors = FALSE
)

# Init the database
shinymanager::create_db(
  credentials_data = credentials,
  sqlite_path = 'credentials.sqlite', # will be created
  passphrase = 'admin_passphrase'
)
```

This will create an sqlite file with an admin user 'shinymanager'.
Next, run carnation pointing to the credentials you just created.

```
run_carnation(credentials='credentials.sqlite', passphrase='admin_passphrase')
```

Now, you will see a login page asking you to authenticate before you can access carnation.
