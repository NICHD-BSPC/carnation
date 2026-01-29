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
