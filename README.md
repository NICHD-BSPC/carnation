# Carnation

**Deeply explore your bulk RNA-Seq data with interactive visualizations**

Carnation is an interactive Shiny dashboard that transforms complex bulk RNA-Seq data into beautiful, insightful visualizations. Designed for both computational and experimental biologists, Carnation makes exploring differential expression analysis, functional enrichment, and pattern analysis intuitive and exciting.

## ✨ Key Features

- **Interactive Visualizations**: Explore your data through multiple perspectives
  - 📊 PCA plots - Visualize sample relationships with gene loadings overlay
  - 🔥 Heatmaps - Examine expression patterns across samples and conditions
  - 📈 MA plots - Identify differentially expressed genes
  - 🧩 UpSet plots - Discover overlapping gene sets across comparisons
  - 📉 Gene plots - Create customizable expression visualizations
- **Functional Enrichment Analysis**: Understand the biological significance of your results
  - 🌐 Network plots - Visualize gene-term relationships
  - 📊 Summary overviews - Compare enrichment across conditions
  - 🔍 Fuzzy search - Find relevant pathways across genes or descriptions
- **Pattern Analysis**: Identify co-regulated gene clusters across conditions
- **Gene Tracking**: Keep your genes of interest at your fingertips with the "Gene scratchpad"
- **Flexible Deployment**: Run locally for personal analysis or on a server to share with collaborators
- **User Management**: Optional authentication system for controlled access in multi-user environments

## 🚀 Installation

### remotes

You can install the developmental version of carnation from github using the `remotes`
package:

```r
install.packages('remotes')
remotes::install_github('NICHD-BSPC/carnation',
                        dependencies=TRUE, build_vignettes=TRUE)
```


### conda

An alternative way to get started with Carnation is through conda, which handles all dependencies automatically:

```bash
# Create environment outside the carnation directory
cd .. && conda env create -p env --file carnation/requirements-pinned.yaml
conda activate ./env
R
```

Then install the package with the `remotes` package. Note, here we set `upgrade='never'`
to make sure the conda-installed package versions remain unchanged.

```r
remotes::install_github('NICHD-BSPC/carnation', upgrade='never')
```

## 🏁 Getting Started

### Data Organization

Organize your data in a directory structure that Carnation can easily navigate:

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

### First Run

Load the library and install required Python dependencies:

```r
library(carnation)
install_carnation()  # Installs plotly and kaleido for PDF export
run_carnation()
```

For remote servers with SSH port forwarding:

```r
run_carnation(options=list(port=12345, launch.browser=FALSE))
```

Then access Carnation at `http://127.0.0.1:12345`

## 🔍 Exploring Your Data

Once your data is loaded, Carnation offers multiple ways to explore:

- **DE Analysis**: Analyze differential expression through multiple visualizations
  - **Summary**: Get a quick overview of your differential expression results
  - **Metadata**: Explore sample metadata and experimental design
  - **PCA Plot**: Visualize sample relationships with gene loadings overlay
  - **MA Plot**: Identify differentially expressed genes with statistical significance
  - **Gene Plot**: Create customizable expression visualizations for genes of interest
  - **UpSet Plot**: Discover overlapping gene sets across multiple comparisons
  - **Heatmap**: Examine expression patterns across samples and conditions

- **Functional Enrichment**: Understand the biological significance of your results
  - **Table**: Interactive tables with powerful search capabilities
  - **Plots**: Seven different visualizations including network plots and dendrograms
  - **Compare Results**: Directly compare enrichment results between conditions

- **Pattern Analysis**: Identify co-regulated gene clusters across conditions
  - **Plot**: Visualize expression patterns of gene clusters
  - **Cluster Membership**: Explore which genes belong to which clusters

- **Gene Scratchpad**: Track genes of interest across all visualizations

## 📚 Documentation

Each module includes comprehensive help documentation accessible through the help buttons throughout the interface. The documentation provides detailed explanations of plot options, statistical methods, and interpretation guidelines.

## 🤝 Contributing

We welcome contributions to Carnation! Please feel free to submit issues or pull requests to the GitHub repository.

## 📄 License

Carnation is available under the MIT license.

## 💻 Server Mode

Carnation supports multi-user environments with authentication:

```r
# Create user database
credentials <- data.frame(
  user = c('shinymanager'),
  password = c('12345'),
  admin = c(TRUE),
  stringsAsFactors = FALSE
)

# Initialize the database
shinymanager::create_db(
  credentials_data = credentials,
  sqlite_path = 'credentials.sqlite',
  passphrase = 'admin_passphrase'
)

# Run with authentication
run_carnation(credentials='credentials.sqlite', passphrase='admin_passphrase')
```
