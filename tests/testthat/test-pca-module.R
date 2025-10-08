library(testthat)
library(shiny)
library(carnation)
library(DESeq2)
library(SummarizedExperiment)
library(plotly)

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
