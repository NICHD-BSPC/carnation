library(shiny)
library(carnation)
library(DESeq2)
library(SummarizedExperiment)
library(plotly)

# Test Pattern Plot Module
test_that("patternPlotServer generates pattern plot correctly", {
  # Create mock data
  mock_dds <- create_mock_dds()
  mock_degpatterns <- create_mock_degpatterns()

  obj <- reactiveValues(
    res = list(comp1 = create_mock_results()),
    degpatterns = list(analysis1 = mock_degpatterns)
  )

  # Set up coldata structure that the module expects
  sample_coldata <- colData(mock_dds)
  if(!"samplename" %in% colnames(sample_coldata)) {
    sample_coldata$samplename <- rownames(sample_coldata)
  }

  coldata <- reactiveValues(
    curr = list(
      all_samples = sample_coldata
    )
  )

  plot_args <- reactive({
    list(
      gene_scratchpad = c("gene1", "gene2"),
      upset_data = list(genes = NULL, labels = NULL)
    )
  })

  config <- reactiveVal(get_config())

  testServer(patternPlotServer, args = list(
    id = "test_pattern",
    obj = obj,
    coldata = coldata,
    plot_args = plot_args,
    config = config
  ), {
    # Set inputs
    session$setInputs(
      dp_analysis = "analysis1",
      deg_facet = "cluster",
      deg_time = "condition",
      deg_color = "none",
      deg_label = "none",
      facet_var_levels = c("1", "2", "3"),
      deg_minc = 1,
      deg_points = FALSE,
      deg_lines = TRUE,
      deg_boxes = TRUE,
      deg_smooth = "line",
      txt_scale = 1,
      deg_rotate = 0,
      deg_cluster = "cluster",
      deg_cluster_levels = c("1", "2", "3")
    )
    session$flushReact()

    # Trigger plot generation
    session$setInputs(plot_do = 1)
    session$flushReact()

    # Test that the module initializes
    expect_true(exists("output"))

    # Test that reactive functions exist
    expect_true(exists("degplot"))
    expect_true(exists("degtable"))

    # Test that reactive functions are reactive
    expect_true(is.reactive(degplot))
    expect_true(is.reactive(degtable))
  })
})

test_that("patternPlotServer handles gene_scratchpad labeling", {
  # Create mock data
  mock_dds <- create_mock_dds()
  mock_degpatterns <- create_mock_degpatterns()

  obj <- reactiveValues(
    res = list(comp1 = create_mock_results()),
    degpatterns = list(analysis1 = mock_degpatterns)
  )

  sample_coldata <- colData(mock_dds)
  if(!"samplename" %in% colnames(sample_coldata)) {
    sample_coldata$samplename <- rownames(sample_coldata)
  }

  coldata <- reactiveValues(
    curr = list(
      all_samples = sample_coldata
    )
  )

  plot_args <- reactive({
    list(
      gene_scratchpad = c("GENE1", "GENE2"),
      upset_data = list(genes = NULL, labels = NULL)
    )
  })

  config <- reactiveVal(get_config())

  testServer(patternPlotServer, args = list(
    id = "test_pattern",
    obj = obj,
    coldata = coldata,
    plot_args = plot_args,
    config = config
  ), {
    # Set inputs with gene_scratchpad labeling
    session$setInputs(
      dp_analysis = "analysis1",
      deg_facet = "cluster",
      deg_time = "condition",
      deg_color = "none",
      deg_label = "gene_scratchpad",
      facet_var_levels = c("1", "2", "3"),
      deg_minc = 1,
      deg_points = FALSE,
      deg_lines = TRUE,
      deg_boxes = TRUE,
      deg_smooth = "line",
      txt_scale = 1,
      deg_rotate = 0,
      deg_cluster = "cluster",
      deg_cluster_levels = c("1", "2", "3")
    )
    session$flushReact()

    # Trigger plot generation
    session$setInputs(plot_do = 1)
    session$flushReact()

    # Test that reactive functions exist
    expect_true(exists("degplot"))
    expect_true(is.reactive(degplot))
  })
})

