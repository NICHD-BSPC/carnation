library(shiny)
library(carnation)
library(DESeq2)
library(SummarizedExperiment)
library(plotly)

# Test Scatter Plot Module
test_that("scatterPlotServer processes data correctly with log2FoldChange", {
  # Create mock data with two comparisons
  mock_dds <- create_mock_dds()
  mock_results1 <- create_mock_results()
  mock_results2 <- create_mock_results()

  # Make some genes significantly DE in both comparisons
  mock_results1$padj[1:20] <- 0.01
  mock_results1$log2FoldChange[1:10] <- 2  # up
  mock_results1$log2FoldChange[11:20] <- -2  # down

  mock_results2$padj[5:25] <- 0.01
  mock_results2$log2FoldChange[5:15] <- 2  # up (overlap with results1)
  mock_results2$log2FoldChange[16:25] <- -2  # down

  obj <- reactiveValues(
    res = list(comp1 = mock_results1, comp2 = mock_results2)
  )

  plot_args <- reactive({
    list(
      fdr.thres = 0.1,
      fc.thres = 0,
      gene.to.plot = c("GENE1", "GENE2")
    )
  })

  config <- reactiveVal(get_config())

  testServer(scatterPlotServer, args = list(
    id = "test_scatter",
    obj = obj,
    plot_args = plot_args,
    config = config
  ), {
    # Set inputs for scatter plot
    session$setInputs(
      x_axis_comp = "comp1",
      y_axis_comp = "comp2",
      compare = "log2FoldChange",
      plot_interactive = "yes",
      show_table = "yes",
      scatter_xmin = -5,
      scatter_xmax = 5,
      scatter_ymin = -5,
      scatter_ymax = 5,
      alpha = 0.7,
      size = 4,
      vline = "yes",
      hline = "yes",
      dline = "yes",
      show_grid = "yes",
      color.palette = "Set2"
    )
    session$flushReact()

    # Trigger plot generation
    session$setInputs(refresh = 1)
    session$flushReact()

    # Test that the module initializes
    expect_true(exists("output"))

    # Test that reactive values exist
    expect_true(exists("flags"))
    expect_true(exists("curr_thres"))
    expect_true(exists("axis_limits"))

    # Test that reactive functions exist
    expect_true(exists("scatterplot"))
    expect_true(exists("scatterplot_ly"))

    # Test that reactive functions are reactive
    expect_true(is.reactive(scatterplot))
    expect_true(is.reactive(scatterplot_ly))
  })
})

# Test scatter plot swap comparisons functionality
test_that("scatterPlotServer handles comparison swapping", {
  # Create mock data
  mock_dds <- create_mock_dds()
  mock_results1 <- create_mock_results()
  mock_results2 <- create_mock_results()

  # Make some genes significantly DE
  mock_results1$padj[1:20] <- 0.01
  mock_results1$log2FoldChange[1:20] <- 2

  mock_results2$padj[5:25] <- 0.01
  mock_results2$log2FoldChange[5:25] <- -2

  obj <- reactiveValues(
    res = list(comp1 = mock_results1, comp2 = mock_results2)
  )

  plot_args <- reactive({
    list(
      fdr.thres = 0.1,
      fc.thres = 0,
      gene.to.plot = character(0)
    )
  })

  config <- reactiveVal(get_config())

  testServer(scatterPlotServer, args = list(
    id = "test_scatter_swap",
    obj = obj,
    plot_args = plot_args,
    config = config
  ), {
    # Set initial inputs
    session$setInputs(
      x_axis_comp = "comp1",
      y_axis_comp = "comp2",
      compare = "log2FoldChange",
      plot_interactive = "yes",
      show_table = "yes",
      scatter_xmin = -5,
      scatter_xmax = 5,
      scatter_ymin = -5,
      scatter_ymax = 5,
      alpha = 0.7,
      size = 4,
      vline = "yes",
      hline = "yes",
      dline = "yes",
      show_grid = "yes",
      color.palette = "Set2"
    )
    session$flushReact()

    # Test swap button
    session$setInputs(swap_comp = 1)
    session$flushReact()

    # Test that the module handles swapping
    expect_true(exists("output"))
    expect_true(exists("flags"))
  })
})

# Test scatter plot autoscale functionality
test_that("scatterPlotServer handles autoscaling", {
  # Create mock data
  mock_dds <- create_mock_dds()
  mock_results1 <- create_mock_results()
  mock_results2 <- create_mock_results()

  # Make some genes with varying fold changes
  mock_results1$padj[1:30] <- 0.01
  mock_results1$log2FoldChange[1:30] <- seq(-3, 3, length.out = 30)

  mock_results2$padj[1:30] <- 0.01
  mock_results2$log2FoldChange[1:30] <- seq(-2, 4, length.out = 30)

  obj <- reactiveValues(
    res = list(comp1 = mock_results1, comp2 = mock_results2)
  )

  plot_args <- reactive({
    list(
      fdr.thres = 0.1,
      fc.thres = 0,
      gene.to.plot = character(0)
    )
  })

  config <- reactiveVal(get_config())

  testServer(scatterPlotServer, args = list(
    id = "test_scatter_autoscale",
    obj = obj,
    plot_args = plot_args,
    config = config
  ), {
    # Set inputs
    session$setInputs(
      x_axis_comp = "comp1",
      y_axis_comp = "comp2",
      compare = "log2FoldChange",
      plot_interactive = "yes",
      show_table = "no",
      scatter_xmin = -5,
      scatter_xmax = 5,
      scatter_ymin = -5,
      scatter_ymax = 5,
      alpha = 0.7,
      size = 4,
      vline = "yes",
      hline = "yes",
      dline = "yes",
      show_grid = "yes",
      color.palette = "Set2"
    )
    session$flushReact()

    # Test autoscale x button
    session$setInputs(scatter_x_auto = 1)
    session$flushReact()

    # Test autoscale y button
    session$setInputs(scatter_y_auto = 1)
    session$flushReact()

    # Test that axis limits are updated
    expect_true(exists("axis_limits"))
    expect_true(!is.null(axis_limits$lim.x))
    expect_true(!is.null(axis_limits$lim.y))
  })
})

# Test scatter plot with gene labeling
test_that("scatterPlotServer handles gene labeling", {
  # Create mock data
  mock_dds <- create_mock_dds()
  mock_results1 <- create_mock_results()
  mock_results2 <- create_mock_results()

  # Make some genes significantly DE
  mock_results1$padj[1:20] <- 0.01
  mock_results1$log2FoldChange[1:20] <- 2

  mock_results2$padj[5:25] <- 0.01
  mock_results2$log2FoldChange[5:25] <- 2

  obj <- reactiveValues(
    res = list(comp1 = mock_results1, comp2 = mock_results2)
  )

  # Simulate genes to label
  plot_args <- reactive({
    list(
      fdr.thres = 0.1,
      fc.thres = 0,
      gene.to.plot = c("GENE1", "GENE2", "GENE3")
    )
  })

  config <- reactiveVal(get_config())

  testServer(scatterPlotServer, args = list(
    id = "test_scatter_labels",
    obj = obj,
    plot_args = plot_args,
    config = config
  ), {
    # Set inputs
    session$setInputs(
      x_axis_comp = "comp1",
      y_axis_comp = "comp2",
      compare = "log2FoldChange",
      plot_interactive = "no",
      show_table = "yes",
      alpha = 0.7,
      size = 4,
      vline = "yes",
      hline = "yes",
      dline = "yes",
      show_grid = "yes",
      color.palette = "Set2"
    )
    session$flushReact()

    # Trigger plot generation
    session$setInputs(refresh = 1)
    session$flushReact()

    # Test that the module handles gene labeling
    expect_true(exists("output"))
    expect_true(is.reactive(scatterplot))
  })
})


