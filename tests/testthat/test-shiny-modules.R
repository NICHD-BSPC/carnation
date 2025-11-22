library(shiny)
library(carnation)
library(DESeq2)
library(SummarizedExperiment)
library(plotly)

# Helper function to create mock DESeq2 data for testing
create_mock_dds <- function(n_genes = 100, n_samples = 6) {
  # Create count matrix
  counts <- matrix(
    rpois(n_genes * n_samples, lambda = 100),
    nrow = n_genes,
    ncol = n_samples
  )
  rownames(counts) <- paste0("gene", 1:n_genes)
  colnames(counts) <- paste0("sample", 1:n_samples)

  # Create sample metadata
  coldata <- data.frame(
    condition = factor(rep(c("control", "treatment"), each = 3)),
    batch = factor(rep(c("A", "B"), times = 3)),
    samplename = colnames(counts),
    row.names = colnames(counts)
  )

  # Create DESeqDataSet
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = coldata,
    design = ~ condition
  )

  return(dds)
}

# Helper function to create mock results
create_mock_results <- function(n_genes = 100) {
  data.frame(
    baseMean = runif(n_genes, 10, 1000),
    log2FoldChange = rnorm(n_genes, 0, 2),
    lfcSE = runif(n_genes, 0.1, 0.5),
    stat = rnorm(n_genes, 0, 3),
    pvalue = runif(n_genes, 0, 1),
    padj = runif(n_genes, 0, 1),
    symbol = paste0("GENE", 1:n_genes),
    row.names = paste0("gene", 1:n_genes)
  )
}

mock_carnation_obj <- function(){
  dds <- create_mock_dds()
  res <- create_mock_results()
  rld <- varianceStabilizingTransformation(dds, blind=TRUE)

  dds_list <- list(main=dds)
  res_list <- list(comp1=res)
  rld_list <- list(main=rld)
  all_dds <- dds
  all_rld <- rld
  dds_mapping <- list(comp1="main")

  list(
    res=res_list,
    dds=dds_list,
    rld=rld_list,
    all_dds=all_dds,
    all_rld=all_rld,
    dds_mapping=dds_mapping,
    enrich=NULL,
    genetonic=NULL,
    degpatterns=NULL
  )
}

# Test PCA Plot Module
test_that("pcaPlotServer handles reactive inputs correctly", {
  # Create mock data
  mock_dds <- create_mock_dds()
  mock_rld <- varianceStabilizingTransformation(mock_dds, blind=TRUE)

  # Create reactive values to simulate app state
  obj <- reactiveValues(
    dds = list(main = mock_dds),
    rld = list(main = mock_rld),
    all_dds = mock_dds,
    all_rld = mock_rld,
    dds_mapping = list(comp1 = 'main')
  )

  # Set up coldata structure that the module expects
  coldata <- reactiveValues(
    curr = list(
      all_samples = colData(mock_dds),
      main = colData(mock_dds)
    )
  )

  config <- reactiveVal(get_config())

  # Test the server function
  testServer(pcaPlotServer, args = list(
    id = "test_pca",
    obj = obj,
    coldata = coldata,
    config = config
  ), {
    # Test that the module initializes correctly
    expect_true(exists("output"))

    # Set inputs in the correct order to avoid reactive dependency issues
    session$setInputs(
      pca_samples = 'all_samples',
      comp_pca = 'main',
      pcx = 1,
      pcy = 2,
      pcz = 'none',
      ntop = 500,
      pca_loadings = 'no',
      pca_loadings_ngenes = 10,
      pca_color = "condition",
      pca_cols = "condition",
      pca_col_levels = c("treatment", "control")
    )
    session$flushReact()

    # Finally trigger the plot
    session$setInputs(
      plot_do = 1
    )

    # Test that reactive values are updated
    expect_true(is.reactive(pca_plot))

    # Test that the module has the expected reactive structure
    expect_true(exists("coldata_pca"))
    expect_true(exists("plot_data"))
    expect_true(exists("data_loaded"))
    expect_is(pca_plot(), 'plotly')
  })
})

# Test the underlying PCA function directly
test_that("plotPCA.san function works correctly", {
  # Create mock data
  mock_dds <- create_mock_dds()
  mock_rld <- varianceStabilizingTransformation(mock_dds, blind=TRUE)

  # Test basic PCA plot generation
  p <- plotPCA.san(mock_rld, intgroup = "condition", pcx = 1, pcy = 2, ntop = 100)

  # Test that a plotly object is returned
  expect_true(heatmaply::is.plotly(p))

  # Test with 3D plot
  p3d <- plotPCA.san(mock_rld, intgroup = "condition", pcx = 1, pcy = 2, pcz = 3, ntop = 100)
  expect_true(heatmaply::is.plotly(p3d))

  # Test with gene loadings
  p_loadings <- plotPCA.san(mock_rld, intgroup = "condition", pcx = 1, pcy = 2,
                           ntop = 100, loadings = TRUE, loadings_ngenes = 5)
  expect_true(heatmaply::is.plotly(p_loadings))
})

# Simplified PCA module test that avoids complex reactive setup
test_that("pcaPlotServer initializes without errors", {
  # Create minimal mock data
  mock_dds <- create_mock_dds()
  mock_rld <- varianceStabilizingTransformation(mock_dds, blind=TRUE)

  # Create minimal reactive values
  obj <- reactiveValues(
    dds = list(main = mock_dds),
    rld = list(main = mock_rld),
    all_dds = mock_dds,
    all_rld = mock_rld,
    dds_mapping = list(comp1 = 'main')
  )

  coldata <- reactiveValues(
    curr = list(
      all_samples = colData(mock_dds),
      main = colData(mock_dds)
    )
  )

  config <- reactiveVal(get_config())

  # Test that the module can be initialized without errors
  testServer(pcaPlotServer, args = list(
    id = "simple_pca",
    obj = obj,
    coldata = coldata,
    config = config
  ), {
    # Just test that the module initializes
    expect_true(exists("output"))
    expect_true(exists("input"))

    # Test that basic reactive values exist
    expect_true(exists("coldata_pca"))
    expect_true(exists("plot_data"))
    expect_true(exists("data_loaded"))
  })
})

## Test MA Plot Module
#test_that("maPlotServer processes data correctly", {
#  # Create mock data
#  mock_results <- create_mock_results()
#
#  obj <- reactiveValues(
#    res = list(test = mock_results)
#  )
#
#  testServer(maPlotServer, args = list(
#    id = "test_ma",
#    obj = reactive(obj),
#    plot_args = reactive(list(
#      fdr.thres = 0.05,
#      fc.thres = 1,
#      gene.to.plot = c("GENE1", "GENE2")
#    ))
#  ), {
#    # Simulate user inputs
#    session$setInputs(
#      de_comp = "test",
#      plot_do = 1
#    )
#
#    # Test that the module responds to inputs
#    expect_true(exists("output"))
#
#    # Test data processing logic
#    # You would add specific assertions based on your module's behavior
#  })
#})
#
## Test Gene Plot Module
#test_that("genePlotServer handles gene selection correctly", {
#  # Create mock data
#  mock_dds <- create_mock_dds()
#
#  obj <- reactiveValues(
#    dds = list(test = mock_dds)
#  )
#
#  testServer(genePlotServer, args = list(
#    id = "test_gene",
#    obj = reactive(obj),
#    coldata = reactive(list(curr = colData(mock_dds))),
#    plot_args = reactive(list(comparison = "test"))
#  ), {
#    # Simulate gene selection
#    session$setInputs(
#      genes_selected = c("gene1", "gene2"),
#      intgroup = "condition",
#      plot_do = 1
#    )
#
#    # Test that gene data is processed correctly
#    expect_true(exists("output"))
#
#    # You could test specific aspects of gene count extraction here
#  })
#})
#
## Test Heatmap Module
#test_that("heatmapServer generates heatmap data correctly", {
#  # Create mock data
#  mock_dds <- create_mock_dds()
#  mock_results <- create_mock_results()
#
#  obj <- reactiveValues(
#    dds = list(test = mock_dds),
#    res = list(test = mock_results)
#  )
#
#  testServer(heatmapServer, args = list(
#    id = "test_heatmap",
#    obj = reactive(obj),
#    coldata = reactive(list(curr = colData(mock_dds))),
#    plot_args = reactive(list(
#      fdr.thres = 0.05,
#      fc.thres = 1,
#      upset_data = list(genes = character(0))
#    )),
#    gene_scratchpad = reactive(character(0))
#  ), {
#    # Simulate user inputs
#    session$setInputs(
#      de_comp = "test",
#      geneset_type = "de",
#      max_gene_num = 50,
#      hmap_scale = "row",
#      plot_do = 1
#    )
#
#    # Test that the module initializes
#    expect_true(exists("output"))
#
#    # Test heatmap gene selection logic
#    # This would involve testing the get_heatmap_genes reactive
#  })
#})
#
## Test Settings Module
#test_that("settingsServer manages user settings correctly", {
#  testServer(settingsServer, args = list(
#    id = "test_settings",
#    details = reactive(list(username = "testuser")),
#    depth = 2,
#    end_offset = 0,
#    assay_fun = function(x) basename(x)
#  ), {
#    # Test that settings are initialized
#    expect_true(exists("output"))
#
#    # Simulate settings changes
#    session$setInputs(
#      project_depth = 3,
#      project_end_offset = 1
#    )
#
#    # Test that settings are updated correctly
#    # You would add specific assertions for your settings logic
#  })
#})
#
## Test Save Object Module
#test_that("saveServer handles object saving correctly", {
#  # Create mock data
#  mock_obj <- list(
#    dds = list(test = create_mock_dds()),
#    res = list(test = create_mock_results())
#  )
#
#  testServer(saveServer, args = list(
#    id = "test_save",
#    original = reactiveValues(obj = mock_obj, path = "/tmp/test.rds"),
#    current = reactiveValues(dds = mock_obj$dds, res = mock_obj$res),
#    coldata = reactive(data.frame()),
#    pattern = reactive(NULL),
#    username = reactive("testuser")
#  ), {
#    # Test save functionality
#    session$setInputs(save_btn = 1)
#
#    # Test that save operations are handled
#    expect_true(exists("output"))
#
#    # You would test the actual save logic here
#  })
#})
#
## Test Load Data Module
#test_that("loadDataServer handles data loading correctly", {
#  testServer(loadDataServer, args = list(
#    id = "test_load",
#    username = reactive("testuser")
#  ), {
#    # Test data loading interface
#    expect_true(exists("output"))
#
#    # Simulate file upload or selection
#    # This would depend on your specific implementation
#
#    # Test that data validation works correctly
#  })
#})
#
## Test helper functions used in modules
#test_that("Module helper functions work correctly", {
#  # Test get_gene_counts function (used in gene plot module)
#  mock_dds <- create_mock_dds()
#
#  # Test with valid inputs
#  result <- get_gene_counts(mock_dds, gene = "gene1", intgroup = "condition")
#  expect_true(is.data.frame(result))
#  expect_true("count" %in% colnames(result))
#  expect_true("gene" %in% colnames(result))
#  expect_true("condition" %in% colnames(result))
#
#  # Test with multiple genes
#  result <- get_gene_counts(mock_dds, gene = c("gene1", "gene2"), intgroup = "condition")
#  expect_equal(length(unique(result$gene)), 2)
#
#  # Test with invalid gene
#  expect_message(
#    get_gene_counts(mock_dds, gene = "nonexistent_gene", intgroup = "condition"),
#    "Skipping"
#  )
#})
#
## Test reactive data flow between modules
#test_that("Module communication works correctly", {
#  # This tests how modules communicate through reactive values
#  # Create a simple test case for gene scratchpad functionality
#
#  gene_scratchpad <- reactiveValues(genes = character(0))
#
#  # Simulate adding genes to scratchpad
#  gene_scratchpad$genes <- c("gene1", "gene2")
#
#  # Test that other modules can access the scratchpad
#  expect_equal(gene_scratchpad$genes, c("gene1", "gene2"))
#
#  # Test scratchpad updates
#  gene_scratchpad$genes <- c(gene_scratchpad$genes, "gene3")
#  expect_equal(length(gene_scratchpad$genes), 3)
#})
