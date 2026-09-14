library(shiny)
library(carnation)
library(DESeq2)
library(SummarizedExperiment)
library(plotly)

test_that("volcanoPlotServer processes data correctly", {
  # Make mock dataset that the module expects
  mock_dds <- create_mock_dds()
  mock_results <- create_mock_results()
  obj <- reactiveValues(
    dds = list(main = mock_dds),
    res = list(test = mock_results)
  )

  config <- reactiveVal(get_config())

  testServer(volcanoPlotServer, args = list(
    id = "test_volcano",
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
      volcano_xmin = -5,
      volcano_xmax = 5,
      volcano_ymin = 0,
      volcano_ymax = 20,
      color_by = "baseMean",
      volcano_alpha = 0.6,
    )

    # Test that the module responds to inputs
    expect_true(exists("output"))

    session$flushReact()

    # Test that reactive values are updated
    expect_true(is.reactive(volcano_plot))
    expect_is(volcano_plot(), "gg")

    expect_true(is.reactive(volcano_plot_ly))
    expect_is(volcano_plot_ly(), "plotly")
  })
})
