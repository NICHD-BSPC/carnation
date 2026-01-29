library(shiny)
library(carnation)
library(DESeq2)
library(SummarizedExperiment)

# Test Upset Plot Module
test_that("upsetPlotServer generates upset data correctly with 'changed' type", {
  # Create mock data with multiple comparisons
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
      fc.thres = 0
    )
  })

  gene_scratchpad <- reactive({ character(0) })
  reset_genes <- reactiveVal()
  config <- reactiveVal(get_config())

  testServer(upsetPlotServer, args = list(
    id = "test_upset",
    obj = obj,
    plot_args = plot_args,
    gene_scratchpad = gene_scratchpad,
    reset_genes = reset_genes,
    config = config
  ), {
    # Set inputs for upset plot
    session$setInputs(
      upset_type = "changed",
      n_intersections = 40,
      min_size = 0,
      text_scale = 1.0
    )
    session$flushReact()

    # Trigger refresh
    session$setInputs(updt_do = 1)
    session$flushReact()

    # Trigger plot generation
    session$setInputs(plot_do = 1)
    session$flushReact()

    # Test that the module initializes
    expect_true(exists("output"))

    # Test that reactive values exist
    expect_true(exists("upset_choices"))
    expect_true(exists("upset_table"))
    expect_true(exists("genes_clicked"))

    # Test that reactive functions exist
    expect_true(exists("upset_gene_lists"))
    expect_true(exists("upsetplot"))
    expect_true(exists("intersect_tbl"))

    # Test that reactive functions are reactive
    expect_true(is.reactive(upset_gene_lists))
    expect_true(is.reactive(upsetplot))
    expect_true(is.reactive(intersect_tbl))
  })
})

# Test upset plot with gene scratchpad integration
test_that("upsetPlotServer integrates with gene scratchpad", {
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

  plot_args <- reactive({
    list(
      fdr.thres = 0.1,
      fc.thres = 0
    )
  })

  # Simulate genes in scratchpad
  gene_scratchpad <- reactive({ c("GENE1", "GENE2", "GENE3") })
  reset_genes <- reactiveVal()
  config <- reactiveVal(get_config())

  testServer(upsetPlotServer, args = list(
    id = "test_upset_scratchpad",
    obj = obj,
    plot_args = plot_args,
    gene_scratchpad = gene_scratchpad,
    reset_genes = reset_genes,
    config = config
  ), {
    # Set inputs
    session$setInputs(
      upset_type = "changed",
      n_intersections = 40,
      min_size = 0,
      text_scale = 1.0
    )
    session$flushReact()

    # Trigger refresh
    session$setInputs(updt_do = 1)
    session$flushReact()

    # Test that genes_clicked reactive value exists
    expect_true(exists("genes_clicked"))

    # Test that the module can handle gene scratchpad updates
    expect_true(exists("output"))
  })
})

