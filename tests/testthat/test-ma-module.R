library(testthat)
library(shiny)
library(carnation)
library(DESeq2)
library(SummarizedExperiment)
library(plotly)

test_that("maPlotServer processes data correctly", {
  # Create mock data
  mock_dds <- create_mock_dds()
  mock_results <- create_mock_results()

  obj <- reactiveValues(
    dds = list(main = mock_dds),
    res = list(test = mock_results)
  )

  # Set up coldata structure that the module expects
  coldata <- reactiveValues(
    curr = list(
      all_samples = colData(mock_dds),
      main = colData(mock_dds)
    )
  )

  config <- reactiveVal(get_config())

  testServer(maPlotServer, args = list(
    id = "test_ma",
    obj = obj,
    plot_args = reactive(list(
      fdr.thres = 0.05,
      fc.thres = 1,
      gene.to.plot = c("GENE1", "GENE2")
    )),
    config = config
  ), {
    # Simulate user inputs
    session$setInputs(
      comp_all = "test",
      ma_ymax = 5,
      ma_ymin = -5
    )

    # Test that the module responds to inputs
    expect_true(exists("output"))

    session$flushReact()

    # Test that reactive values are updated
    expect_true(is.reactive(maplot))
    expect_is(maplot(), "gg")

    expect_true(is.reactive(maplot_ly))
    expect_is(maplot_ly(), "plotly")
    # Test data processing logic
    # You would add specific assertions based on your module's behavior
  })
})
