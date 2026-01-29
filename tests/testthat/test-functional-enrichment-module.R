library(shiny)
library(carnation)
library(DESeq2)
library(ggplot2)

# Test Functional Enrichment Module
test_that("enrichServer initializes correctly with enrichment data", {
  # Create mock data
  mock_results <- create_mock_results()
  mock_enrich <- create_mock_enrichment()
  mock_genetonic <- create_mock_genetonic(mock_enrich, mock_results)

  obj <- reactiveValues(
    res = list(comp1 = mock_results),
    enrich = list(
      comp1 = list(
        res = "comp1",
        changed = list(
          BP = mock_enrich
        )
      )
    ),
    genetonic = list(
      comp1 = list(
        res = "comp1",
        changed = list(
          BP = mock_genetonic
        )
      )
    )
  )

  upset_table <- reactiveValues(
    tbl = NULL,
    intersections = NULL,
    set_labels = NULL
  )

  gene_scratchpad <- reactive({ c("gene1", "gene2") })
  reset_genes <- reactive({ FALSE })
  config <- reactiveVal(get_config())

  testServer(enrichServer, args = list(
    id = "test_enrich",
    obj = obj,
    upset_table = upset_table,
    gene_scratchpad = gene_scratchpad,
    reset_genes = reset_genes,
    config = config
  ), {
    # Set inputs
    session$setInputs(
      comp_tbl = "comp1",
      geneset_tbl = "changed",
      pathway_tbl = "BP",
      comp_fun = "comp1",
      geneset = "changed",
      pathway = "BP",
      search_opts = "genes + description",
      search_txt = "",
      subset_opts = "none",
      genes.per.line = 6
    )
    session$flushReact()

    # Test that the module initializes
    expect_true(exists("output"))

    # Test that reactive values exist
    expect_true(exists("enrich_data"))
    expect_true(exists("flags"))
    expect_true(exists("genes_clicked"))

    # Test that reactive functions exist
    expect_true(exists("get_func_table"))
    expect_true(is.reactive(get_func_table))

    # Test that distill & fuzzy tables exist
    expect_true(!is.null(enrich_data$distill))
    expect_true(!is.null(enrich_data$fuzzy))

    # Test enrich/genetonic obj (used for plots) exist
    expect_true(is.reactive(enrich_obj))
    expect_true(is.reactive(genetonic_obj))
  })
})

