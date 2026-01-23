library(shiny)
library(carnation)
library(DESeq2)

# Test Load New Data Module
test_that("loadDataServer initializes correctly", {
  username <- "testuser"
  config <- reactiveVal(get_config())

  testServer(loadDataServer, args = list(
    id = "test_load",
    username = username,
    config = config,
    rds = NULL
  ), {
    # Test that the module initializes
    expect_true(exists("output"))

    # Test that reactive values exist
    expect_true(exists("new_obj"))
    expect_true(exists("reload_parent"))

    # Test that new_obj has expected structure
    expect_true(is.null(new_obj$res_list))
    expect_true(is.null(new_obj$dds_list))
    expect_true(is.null(new_obj$rld_list))
    expect_true(is.null(new_obj$enrich_list))
    expect_true(is.null(new_obj$degpatterns))
    expect_true(is.null(new_obj$genetonic))

    # Test that reload_parent is reactive
    expect_true(is.reactive(reload_parent))
    expect_false(reload_parent())
  })
})

test_that("loadDataServer handles existing object for editing", {
  username <- "testuser"
  config <- reactiveVal(get_config())

  # Create mock carnation object
  obj <- mock_carnation_obj()

  rds <- reactiveValues(obj = obj)

  testServer(loadDataServer, args = list(
    id = "test_load",
    username = username,
    config = config,
    rds = rds
  ), {
    # Test that edit_obj reactive exists
    expect_true(exists("edit_obj"))
    expect_true(is.reactive(edit_obj))

    session$flushReact()

    # Test that new_obj is populated from existing object
    expect_true(!is.null(new_obj$res_list))
    expect_true(!is.null(new_obj$dds_list))
    expect_true(!is.null(new_obj$rld_list))
  })
})

