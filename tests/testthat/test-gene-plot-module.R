library(shiny)
library(carnation)
library(DESeq2)
library(SummarizedExperiment)
library(plotly)

# Test Gene Plot Module
test_that("genePlotServer handles gene selection correctly", {
  # Create mock data
  mock_dds <- create_mock_dds()
  mock_rld <- varianceStabilizingTransformation(mock_dds, blind=TRUE)
  mock_results <- create_mock_results()

  # Create reactive values to simulate app state
  obj <- reactiveValues(
    dds = list(main = mock_dds),
    rld = list(main = mock_rld),
    res = list(comp1 = mock_results),
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

  testServer(genePlotServer, args = list(
    id = "test_gene",
    obj = obj,
    coldata = coldata,
    plot_args = reactive(
      list(
        gene.to.plot = c("gene1", "gene2"),
        gene.id = rownames(mock_dds),
        comp_all = "comp1"  # Fixed: should match the key in res list
      )
    ),
    config = config
  ), {
    # Set inputs in the correct order to avoid validation errors
    # First set the basic sample selection
    session$setInputs(
      samples = "all_samples",  # This should match a key in coldata$curr
      xvar = "condition",  # This should match a column in the coldata
      norm_method = "vst",
      ymax = 10,
      ymin = -10,
      logy = FALSE,
      boxes = TRUE,
      freey = TRUE,
      x_rotate = 30,
      color = "batch",
      trendline = "line",
      gene_nrow = 1,
      legend = TRUE,
      txt_scale = 1,
      facet = "none",           # Set to 'none' to avoid faceting issues
      facet_vars = "",          # Empty string to avoid NULL issues
      facet_var_levels = character(0)  # Empty character vector
    )
    session$flushReact()

    # Test that gene data is processed correctly
    expect_true(exists("output"))

    # Test that basic reactive values exist
    expect_true(exists("gene_coldata"))
    expect_true(exists("gene_plot_data"))

    expect_true(any(plot_args()$gene.to.plot %in% plot_args()$gene.id))
    expect_true(any(plot_args()$gene.to.plot %in% rownames(gene_plot_data$all)))

    # Test that the reactive functions exist
    expect_true(is.reactive(normplot))

    expect_true(!is.null(gene_plot_data$plotted))
    expect_true(!is.null(gene_plot_data$handle))

    expect_is(gene_plot_data$handle, "ggplot")
  })
})
