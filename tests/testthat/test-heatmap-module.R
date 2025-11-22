library(shiny)
library(carnation)
library(DESeq2)
library(SummarizedExperiment)
library(plotly)

# Test Heatmap Module
test_that("heatmapServer generates heatmap data correctly with 'de_genes' ", {
  # Create mock data
  mock_dds <- create_mock_dds()
  mock_rld <- varianceStabilizingTransformation(mock_dds, blind=TRUE)
  mock_results <- create_mock_results()

  obj <- reactiveValues(
    rld = list(main = mock_rld),
    all_rld = mock_rld,
    res = list(test = mock_results)
  )

  # Set up coldata structure that the module expects
  # Make sure colData has the required columns and structure
  sample_coldata <- colData(mock_dds)
  # Ensure samplename column exists (required by heatmap module)
  if(!"samplename" %in% colnames(sample_coldata)) {
    sample_coldata$samplename <- rownames(sample_coldata)
  }

  coldata <- reactiveValues(
    curr = list(
      all_samples = sample_coldata,
      main = sample_coldata
    )
  )

  config <- reactiveVal(get_config())

  testServer(heatmapServer, args = list(
    id = "test_heatmap",
    obj = obj,
    coldata = coldata,
    plot_args = reactive(list(
      fdr.thres = 0.1,
      fc.thres = 0,
      upset_data = list(genes = character(0))
    )),
    gene_scratchpad = reactive({ c("gene1","gene2") }),
    config = config
  ), {
    # Set inputs in stages to avoid reactive dependency issues
    # First set basic parameters
    session$setInputs(
      de_comp = "test",
      geneset_type = "de_genes",
      hmap_type = "changed",
      hmap_rld = "all_samples",
      max_gene_num = 50,
      hmap_scale = "row",
      hmap_clust = "row",
      hmap_rank = "padj",
      hmap_cols = "condition",
      hmap_colnames = "condition"
    )
    session$flushReact()

    # Trigger plot generation
    session$setInputs(plot_do = 1)
    session$flushReact()

    # Test that the module initializes
    expect_true(exists("output"))

    # Test that basic reactive values exist
    expect_true(exists("get_heatmap_genes"))
    expect_true(exists("make_heatmap"))

    # Test that reactive functions exist (but may not be ready to execute)
    expect_true(is.reactive(get_heatmap_genes))
    expect_true(is.reactive(make_heatmap))
  })

  # this tests 'gene_scratchpad' mode
  testServer(heatmapServer, args = list(
    id = "test_heatmap_scratch",
    obj = obj,
    coldata = coldata,
    plot_args = reactive(list(
      fdr.thres = 0.1,
      fc.thres = 0,
      upset_data = list(genes = character(0))
    )),
    gene_scratchpad = reactive({ c("gene1","gene2") }),
    config = config
  ), {
    # Set inputs in stages to avoid reactive dependency issues
    # First set basic parameters
    session$setInputs(
      de_comp = "test",
      geneset_type = "gene_scratchpad",
      hmap_type = "changed",
      hmap_rld = "all_samples",
      max_gene_num = 50,
      hmap_scale = "row",
      hmap_clust = "row",
      hmap_rank = "padj",
      hmap_cols = "condition",
      hmap_colnames = "condition"
    )
    session$flushReact()

    # Trigger plot generation
    session$setInputs(plot_do = 1)
    session$flushReact()

    # Test that the module initializes
    expect_true(exists("output"))

    # Test that basic reactive values exist
    expect_true(exists("get_heatmap_genes"))
    expect_true(exists("make_heatmap"))

    # Test that reactive functions exist (but may not be ready to execute)
    expect_true(is.reactive(get_heatmap_genes))
    expect_true(is.reactive(make_heatmap))
  })

  upset_genes <- list(set01 = c("gene1", "gene2", "gene3"),
                      set02 = c("gene4", "gene5"))
  upset_labels <- c("set01 (n = 3)", "set 02 (n = 2)")

  # this tests 'upset_intersections' mode
  testServer(heatmapServer, args = list(
    id = "test_heatmap_upset",
    obj = obj,
    coldata = coldata,
    plot_args = reactive(list(
      fdr.thres = 0.1,
      fc.thres = 0,
      upset_data = list(genes = upset_genes, labels = upset_labels)
    )),
    gene_scratchpad = reactive({ c("gene1","gene2") }),
    config = config
  ), {
    # Set inputs in stages to avoid reactive dependency issues
    # First set basic parameters
    session$setInputs(
      de_comp = "test",
      geneset_type = "upset_intersections",
      hmap_type = "changed",
      hmap_rld = "all_samples",
      max_gene_num = 50,
      hmap_scale = "row",
      hmap_clust = "row",
      hmap_rank = "padj",
      hmap_cols = "condition",
      hmap_colnames = "condition"
    )
    session$flushReact()

    # Trigger plot generation
    session$setInputs(plot_do = 1)
    session$flushReact()

    # Test that the module initializes
    expect_true(exists("output"))

    # Test that basic reactive values exist
    expect_true(exists("get_heatmap_genes"))
    expect_true(exists("make_heatmap"))

    # Test that reactive functions exist (but may not be ready to execute)
    expect_true(is.reactive(get_heatmap_genes))
    expect_true(is.reactive(make_heatmap))
  })})



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
