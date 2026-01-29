library(shiny)
library(carnation)
library(DESeq2)

# Test Save Object Module
test_that("saveServer initializes correctly", {
  # Create mock carnation object
  obj <- mock_carnation_obj()

  # Create temporary path
  temp_dir <- tempdir()
  temp_path <- file.path(temp_dir, "test.rnaseq.rds")

  original <- reactiveValues(
    obj = obj,
    path = temp_path
  )

  current <- reactiveValues(
    res = obj$res,
    dds = obj$dds,
    rld = obj$rld
  )

  coldata <- reactive({
    lapply(obj$dds, colData)
  })

  pattern <- "rnaseq"
  username <- reactive({ "testuser" })
  config <- reactiveVal(get_config())

  testServer(saveServer, args = list(
    id = "test_save",
    original = original,
    current = current,
    coldata = coldata,
    pattern = pattern,
    username = username,
    config = config
  ), {
    # Test that the module initializes
    expect_true(exists("output"))

    # Test that reactive values exist
    expect_true(exists("reload_parent"))
    expect_true(exists("save_flag"))

    # Test initial state
    expect_false(reload_parent$flag)
    expect_false(save_flag$l)
  })
})

test_that("saveServer handles NULL original object", {
  original <- reactiveValues(
    obj = NULL,
    path = NULL
  )

  current <- reactiveValues()
  coldata <- reactive({ list() })
  pattern <- "rnaseq"
  username <- reactive({ "testuser" })
  config <- reactiveVal(get_config())

  testServer(saveServer, args = list(
    id = "test_save",
    original = original,
    current = current,
    coldata = coldata,
    pattern = pattern,
    username = username,
    config = config
  ), {
    # Test that reactive values exist even with NULL object
    expect_true(exists("reload_parent"))
    expect_true(exists("save_flag"))

    # Test initial state
    expect_false(reload_parent$flag)
    expect_false(save_flag$l)
  })
})

