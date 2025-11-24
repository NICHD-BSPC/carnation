library(carnation)
library(dplyr)
library(shiny)
library(DESeq2)
library(SummarizedExperiment)
library(ggplot2)
library(ggrepel)
library(plotly)

# Test get_y_init function
test_that("get_y_init calculates correct y-axis limits", {
  # Create test data frame
  df <- data.frame(count = c(10, 20, 30, 40, 50))

  # Test with y_delta = 0.1 and pseudocount = 1
  result <- get_y_init(df, y_delta = 0.1, pseudocount = 1)

  # Expected values:
  # df$count becomes c(11, 21, 31, 41, 51)
  # df.max = 51, df.min = 11
  # max.init = round((1 + 0.1) * 51) = round(56.1) = 56
  # min.init = round((1 - 0.1) * 11) = round(9.9) = 10
  expect_equal(result, c(10, 56))

  # Test with y_delta = 0.2 and pseudocount = 0
  result <- get_y_init(df, y_delta = 0.2, pseudocount = 0)

  # Expected values:
  # df$count remains c(10, 20, 30, 40, 50)
  # df.max = 50, df.min = 10
  # max.init = round((1 + 0.2) * 50) = round(60) = 60
  # min.init = round((1 - 0.2) * 10) = round(8) = 8
  expect_equal(result, c(8, 60))

  # Test with very small values and pseudocount
  df_small <- data.frame(count = c(0.1, 0.2, 0.3))
  result <- get_y_init(df_small, y_delta = 0.1, pseudocount = 0.5)

  # Expected values:
  # df$count becomes c(0.6, 0.7, 0.8)
  # df.max = 0.8, df.min = 0.6
  # max.init = round((1 + 0.1) * 0.8) = round(0.88) = 1
  # min.init = round((1 - 0.1) * 0.6) = round(0.54) = 1
  # But since min.init < 1, it becomes pseudocount*0.1 = 0.5*0.1 = 0.05
  expect_equal(result, c(0.05, 1))
})

test_that("get_y_init throws error for invalid inputs", {
  # Test with missing 'count' column
  df_invalid <- data.frame(values = c(10, 20, 30))
  expect_error(get_y_init(df_invalid, y_delta = 0.1, pseudocount = 1),
               "Column \"count\" not found in data frame")

  # Test with invalid y_delta (negative)
  df_valid <- data.frame(count = c(10, 20, 30))
  expect_error(get_y_init(df_valid, y_delta = -0.1, pseudocount = 1),
               "y_delta must be between 0 and 1")

  # Test with invalid y_delta (greater than 1)
  expect_error(get_y_init(df_valid, y_delta = 1.1, pseudocount = 1),
               "y_delta must be between 0 and 1")
})

# Test add.set.column function
test_that("add.set.column adds set column correctly", {
  # Create a test matrix
  test_df <- data.frame(
    symbol = c("gene1", "gene2", "gene3"),
    set1 = c(1, 0, 0),
    set2 = c(1, 1, 0),
    set3 = c(0, 1, 1)
  )

  result <- add.set.column(test_df)

  # Check that the result has the correct columns
  expect_true("set" %in% colnames(result))
  expect_true("symbol" %in% colnames(result))

  # Check that the set column is at the beginning
  expect_equal(colnames(result)[2], "set")
  expect_equal(colnames(result)[1], "symbol")

  # Check that the set values are assigned correctly
  # The function should group rows with identical patterns
  # and assign them set names like "set01", "set02", etc.
  unique_sets <- unique(result$set)
  expect_true(length(unique_sets) <= 3) # At most 3 unique patterns in our test data

  # Check that rows with the same pattern have the same set value
  pattern1 <- paste(c(1, 1, 0), collapse = " ")
  pattern2 <- paste(c(0, 1, 1), collapse = " ")
  pattern3 <- paste(c(0, 0, 1), collapse = " ")

  # Get rows with each pattern
  rows_pattern1 <- which(apply(test_df[, c("set1", "set2", "set3")], 1, function(x) paste(x, collapse = " ")) == pattern1)
  if (length(rows_pattern1) > 0) {
    set_value <- result$set[rows_pattern1[1]]
    expect_true(all(result$set[rows_pattern1] == set_value))
  }
})

# Test format_genes function
test_that("format_genes formats gene lists correctly", {
  # Test with a single gene
  expect_equal(format_genes("gene1", sep="\\/", genes.per.line=6), "gene1")

  # Test with multiple genes but fewer than genes.per.line
  gene_list <- "gene1/gene2/gene3"
  expected <- "gene1,gene2,gene3"
  expect_equal(format_genes(gene_list, sep="\\/", genes.per.line=6), expected)

  # Test with more genes than genes.per.line
  gene_list <- "gene1/gene2/gene3/gene4/gene5/gene6/gene7/gene8"
  # Expected: first 6 genes in one group, then the rest
  expected <- "gene1,gene2,gene3,gene4,gene5,gene6, gene7,gene8,"
  expect_equal(format_genes(gene_list, sep="\\/", genes.per.line=6), expected)

  # Test with multiple gene lists
  gene_lists <- c(
    "gene1/gene2/gene3",
    "gene4/gene5/gene6/gene7/gene8/gene9/gene10"
  )
  expected <- c(
    "gene1,gene2,gene3",
    "gene4,gene5,gene6,gene7,gene8,gene9, gene10,"
  )
  expect_equal(format_genes(gene_lists, sep="\\/", genes.per.line=6), expected)
})

# Test top.genes function
test_that("top.genes returns correct top genes", {
  # Create a mock DESeq2 results data frame
  res <- data.frame(
    baseMean = 1:10,
    log2FoldChange = c(3, 2.5, 2, 1.5, 1, -1, -1.5, -2, -1.5, -2),
    padj = c(0.001, 0.005, 0.01, 0.02, 0.05, 0.05, 0.02, 0.01, 0.005, 0.001),
    symbol = paste0("gene", 1:10),
    row.names = paste0("ENSG", 1:10)
  )

  # Test with default parameters (by log2FoldChange)
  result <- top.genes(res, fdr.thres=0.05, fc.thres=0, n=4)
  expect_equal(length(result), 4)
  expect_true(all(result %in% paste0("gene", c(1, 2, 8, 10))))

  # Test with by='padj'
  result <- top.genes(res, fdr.thres=0.05, fc.thres=0, n=4, by='padj')
  expect_equal(length(result), 4)
  expect_true(all(result %in% paste0("gene", c(1, 2, 9, 10))))

  # Test with fc.thres > 0
  result <- top.genes(res, fdr.thres=0.05, fc.thres=2, n=4)
  expect_equal(length(result), 4)
  expect_true(all(result %in% paste0("gene", c(1, 2, 8, 10))))

  # Test with n greater than number of DE genes
  result <- top.genes(res, fdr.thres=0.005, fc.thres=0, n=4)
  expect_equal(length(result), 2)
  expect_true(all(result %in% paste0("gene", c(1, 10))))
})

# Test get_project_name_from_path function
test_that("get_project_name_from_path extracts project name correctly", {
  # Test with standard path
  path <- "/path/to/project/test/main.pattern.rds"
  result <- get_project_name_from_path(path, depth=2, end_offset=0, fsep="/")
  expect_equal(result, "project/test")

  # Test with staging directory
  path <- "/path/to/project/dev/test/main.pattern.rds"
  result <- get_project_name_from_path(path, depth=2, end_offset=0, staging_dir="dev", fsep="/")
  expect_equal(result, "project/dev/test")

  # Test with end_offset
  path <- "/path/to/project/test/subdir/main.pattern.rds"
  result <- get_project_name_from_path(path, depth=3, end_offset=1, fsep="/")
  expect_equal(result, "project/test")
})

# Test fromList.with.names function
test_that("fromList.with.names creates correct binary matrix", {
  # Create test list
  test_list <- list(
    set1 = c(gene1 = "Gene1", gene2 = "Gene2", gene3 = "Gene3"),
    set2 = c(gene2 = "Gene2", gene3 = "Gene3", gene4 = "Gene4"),
    set3 = c(gene3 = "Gene3", gene4 = "Gene4", gene5 = "Gene5")
  )

  result <- fromList.with.names(test_list)

  # Check dimensions and column names
  expect_equal(ncol(result), 4)  # 3 sets + symbol column
  expect_equal(nrow(result), 5)  # 5 unique genes
  expect_equal(colnames(result)[2:4], c("set1", "set2", "set3"))

  # Check that the first column is 'symbol'
  expect_equal(colnames(result)[1], "symbol")

  # Check binary values for specific genes
  # gene1 should be in set1 only
  gene1_row <- which(rownames(result) == "gene1")
  expect_equal(as.numeric(result[gene1_row, c("set1", "set2", "set3")]), c(1, 0, 0))

  # gene3 should be in all sets
  gene3_row <- which(rownames(result) == "gene3")
  expect_equal(as.numeric(result[gene3_row, c("set1", "set2", "set3")]), c(1, 1, 1))

  # gene5 should be in set3 only
  gene5_row <- which(rownames(result) == "gene5")
  expect_equal(as.numeric(result[gene5_row, c("set1", "set2", "set3")]), c(0, 0, 1))

  # Check symbol values
  expect_equal(result$symbol[gene1_row], "Gene1")
  expect_equal(result$symbol[gene3_row], "Gene3")
})

# Test PCA plotting function
test_that("plotPCA.san function works correctly", {
  # Create mock data
  mock_dds <- create_mock_dds()
  mock_rld <- varianceStabilizingTransformation(mock_dds, blind=TRUE)

  # Test basic PCA plot generation
  p <- plotPCA.san(mock_rld, intgroup = "condition", pcx = 1, pcy = 2, ntop = 100)

  # Test that a plotly object is returned
  expect_true(heatmaply::is.plotly(p))

  # Test with 3D plot
  p3d <- plotPCA.san(mock_rld, intgroup = "condition", pcx = 1, pcy = 2, pcz = 3, ntop = 100)
  expect_true(heatmaply::is.plotly(p3d))

  # Test with gene loadings
  p_loadings <- plotPCA.san(mock_rld, intgroup = "condition", pcx = 1, pcy = 2,
                           ntop = 100, loadings = TRUE, loadings_ngenes = 5)
  expect_true(heatmaply::is.plotly(p_loadings))
})

# Test plotMA.label function
test_that("plotMA.label creates correct MA plot", {
  # Create mock results data
  mock_res <- create_mock_results(50)

  # Make some genes significantly differentially expressed
  mock_res$padj[1:10] <- 0.001  # Significant genes
  mock_res$log2FoldChange[1:5] <- c(2, 3, -2, -3, 2.5)  # Up and down regulated
  mock_res$log2FoldChange[6:10] <- c(-1.5, 1.8, -2.2, 2.8, -1.2)

  # Test basic MA plot generation
  p <- plotMA.label(mock_res, fdr.thres = 0.01, fc.thres = 0)

  # Test that a ggplot object is returned
  expect_true(is_ggplot(p))

  # Test that the plot has the expected layers
  expect_true(length(p$layers) > 0)

  # Test with custom parameters
  p_custom <- plotMA.label(mock_res, fdr.thres = 0.05, fc.thres = 1,
                          fc.lim = c(-4, 4))
  expect_true(is_ggplot(p_custom))

  # Test with gene labels
  lab_genes <- c("GENE1", "GENE2", "GENE3")
  p_labeled <- plotMA.label(mock_res, fdr.thres = 0.01, fc.thres = 0,
                           lab.genes = lab_genes)
  expect_true(is_ggplot(p_labeled))

  # Test that the plot has more layers when genes are labeled
  expect_true(length(p_labeled$layers) >= length(p$layers))
})

test_that("plotMA.label handles edge cases correctly", {
  # Test with minimal data
  minimal_res <- data.frame(
    baseMean = c(100, 200, 300),
    log2FoldChange = c(1, -1, 0),
    padj = c(0.001, 0.05, 0.5),
    symbol = c("GENE1", "GENE2", "GENE3"),
    row.names = c("gene1", "gene2", "gene3")
  )

  p <- plotMA.label(minimal_res)
  expect_true(is_ggplot(p))

  # Test with missing symbol column
  no_symbol_res <- minimal_res
  no_symbol_res$symbol <- NULL

  p_no_symbol <- plotMA.label(no_symbol_res)
  expect_true(is_ggplot(p_no_symbol))

  # Test with NA values
  na_res <- minimal_res
  na_res$log2FoldChange[1] <- NA
  na_res$padj[2] <- NA

  p_na <- plotMA.label(na_res)
  expect_true(is_ggplot(p_na))
})

test_that("plotMA.label throws errors for invalid inputs", {
  # Test with missing required columns
  invalid_res <- data.frame(
    baseMean = c(100, 200),
    pvalue = c(0.01, 0.05)
  )

  expect_error(plotMA.label(invalid_res),
               'DE analysis results must contain "padj" & "log2FoldChange" columns')

  # Test with missing padj column
  missing_padj <- data.frame(
    baseMean = c(100, 200),
    log2FoldChange = c(1, -1)
  )

  expect_error(plotMA.label(missing_padj),
               'DE analysis results must contain "padj" & "log2FoldChange" columns')
})

# Test plotMA.label_ly function (interactive plotly version)
test_that("plotMA.label_ly creates correct interactive MA plot", {
  # Create mock results data
  mock_res <- create_mock_results(50)

  # Make some genes significantly differentially expressed
  mock_res$padj[1:10] <- 0.001  # Significant genes
  mock_res$log2FoldChange[1:5] <- c(2, 3, -2, -3, 2.5)  # Up and down regulated
  mock_res$log2FoldChange[6:10] <- c(-1.5, 1.8, -2.2, 2.8, -1.2)

  # Test basic interactive MA plot generation
  p <- plotMA.label_ly(mock_res, fdr.thres = 0.01, fc.thres = 0)

  # Test that a plotly object is returned
  expect_true(heatmaply::is.plotly(p))

  # Test with custom parameters
  p_custom <- plotMA.label_ly(mock_res, fdr.thres = 0.05, fc.thres = 1,
                             fc.lim = c(-4, 4))
  expect_true(heatmaply::is.plotly(p_custom))

  # Test with gene labels
  lab_genes <- c("GENE1", "GENE2", "GENE3")
  p_labeled <- plotMA.label_ly(mock_res, fdr.thres = 0.01, fc.thres = 0,
                              lab.genes = lab_genes)
  expect_true(heatmaply::is.plotly(p_labeled))
})

test_that("plotMA.label_ly handles edge cases correctly", {
  # Test with minimal data
  minimal_res <- data.frame(
    baseMean = c(100, 200, 300),
    log2FoldChange = c(1, -1, 0),
    padj = c(0.001, 0.05, 0.5),
    symbol = c("GENE1", "GENE2", "GENE3"),
    row.names = c("gene1", "gene2", "gene3")
  )

  p <- plotMA.label_ly(minimal_res)
  expect_true(heatmaply::is.plotly(p))

  # Test with missing symbol column
  no_symbol_res <- minimal_res
  no_symbol_res$symbol <- NULL

  p_no_symbol <- plotMA.label_ly(no_symbol_res)
  expect_true(heatmaply::is.plotly(p_no_symbol))

  # Test with NA values
  na_res <- minimal_res
  na_res$log2FoldChange[1] <- NA
  na_res$padj[2] <- NA

  p_na <- plotMA.label_ly(na_res)
  expect_true(heatmaply::is.plotly(p_na))
})

test_that("plotMA.label_ly throws errors for invalid inputs", {
  # Test with missing required columns
  invalid_res <- data.frame(
    baseMean = c(100, 200),
    pvalue = c(0.01, 0.05)
  )

  expect_error(plotMA.label_ly(invalid_res),
               'DE analysis results must contain "padj" & "log2FoldChange" columns')

  # Test with missing log2FoldChange column
  missing_lfc <- data.frame(
    baseMean = c(100, 200),
    padj = c(0.01, 0.05)
  )

  expect_error(plotMA.label_ly(missing_lfc),
               'DE analysis results must contain "padj" & "log2FoldChange" columns')
})

# Test MA plot data processing logic
test_that("MA plot functions process data correctly", {
  # Create test data with known values
  test_res <- data.frame(
    baseMean = c(100, 200, 300, 400, 500),
    log2FoldChange = c(2, -2, 0.5, -0.5, 10),  # Include extreme value
    padj = c(0.001, 0.001, 0.1, 0.1, 0.001),
    symbol = c("UP1", "DOWN1", "NS1", "NS2", "EXTREME"),
    row.names = c("gene1", "gene2", "gene3", "gene4", "gene5")
  )

  # Test that extreme values are handled (should be clipped to fc.lim)
  p_ggplot <- plotMA.label(test_res, fdr.thres = 0.01, fc.thres = 0, fc.lim = c(-5, 5))
  expect_true(is_ggplot(p_ggplot))

  p_plotly <- plotMA.label_ly(test_res, fdr.thres = 0.01, fc.thres = 0, fc.lim = c(-5, 5))
  expect_true(heatmaply::is.plotly(p_plotly))

  # Test with different significance thresholds
  p_strict <- plotMA.label(test_res, fdr.thres = 0.0001, fc.thres = 1)
  expect_true(is_ggplot(p_strict))

  p_lenient <- plotMA.label(test_res, fdr.thres = 0.1, fc.thres = 0)
  expect_true(is_ggplot(p_lenient))
})

# Test symbol column handling
test_that("MA plot functions handle symbol columns correctly", {
  # Test with SYMBOL column (uppercase)
  res_symbol <- data.frame(
    baseMean = c(100, 200),
    log2FoldChange = c(1, -1),
    padj = c(0.01, 0.05),
    SYMBOL = c("GENE1", "GENE2"),
    row.names = c("gene1", "gene2")
  )

  p1 <- plotMA.label(res_symbol)
  expect_true(is_ggplot(p1))

  p1_ly <- plotMA.label_ly(res_symbol)
  expect_true(heatmaply::is.plotly(p1_ly))

  # Test with ALIAS column
  res_alias <- data.frame(
    baseMean = c(100, 200),
    log2FoldChange = c(1, -1),
    padj = c(0.01, 0.05),
    ALIAS = c("ALIAS1", "ALIAS2"),
    row.names = c("gene1", "gene2")
  )

  p2 <- plotMA.label(res_alias)
  expect_true(is_ggplot(p2))

  p2_ly <- plotMA.label_ly(res_alias)
  expect_true(heatmaply::is.plotly(p2_ly))

  # Test with NA symbols
  res_na_symbol <- data.frame(
    baseMean = c(100, 200),
    log2FoldChange = c(1, -1),
    padj = c(0.01, 0.05),
    symbol = c("GENE1", NA),
    row.names = c("gene1", "gene2")
  )

  p3 <- plotMA.label(res_na_symbol)
  expect_true(is_ggplot(p3))

  p3_ly <- plotMA.label_ly(res_na_symbol)
  expect_true(heatmaply::is.plotly(p3_ly))
})

# Helper function to create mock gene count data for getcountplot
create_mock_gene_counts <- function(n_samples = 6, genes = c("gene1", "gene2")) {
  # Create sample metadata
  conditions <- factor(rep(c("control", "treatment"), each = n_samples/2))
  batches <- factor(rep(c("A", "B"), times = n_samples/2))
  samples <- paste0("sample", 1:n_samples)

  # Create count data for each gene
  count_data <- lapply(genes, function(gene) {
    data.frame(
      count = rpois(n_samples, lambda = 100) + runif(n_samples, 50, 150),
      gene = gene,
      condition = conditions,
      batch = batches,
      sample = samples,
      stringsAsFactors = FALSE
    )
  })

  # Combine all genes into one data frame
  do.call(rbind, count_data)
}

# Test getcountplot function
test_that("getcountplot creates correct gene count plots", {
  # Create mock gene count data
  mock_counts <- create_mock_gene_counts(n_samples = 6, genes = c("gene1", "gene2"))

  # Test basic plot generation
  p <- getcountplot(
    df = mock_counts,
    intgroup = "condition",
    factor.levels = c("control", "treatment"),
    color = "gene"
  )

  # Test that a ggplot object is returned
  expect_true(is.ggplot(p))

  # Test that the plot has the expected layers
  expect_true(length(p$layers) > 0)

  # Test with custom parameters
  p_custom <- getcountplot(
    df = mock_counts,
    intgroup = "condition",
    factor.levels = c("control", "treatment"),
    color = "gene",
    ymin = 0,
    ymax = 500,
    log = FALSE,
    boxes = FALSE,
    title = "Test Gene Plot"
  )

  expect_true(is.ggplot(p_custom))
  expect_equal(p_custom$labels$title, "Test Gene Plot")
})

test_that("getcountplot handles different plot options correctly", {
  # Create mock data
  mock_counts <- create_mock_gene_counts(n_samples = 8, genes = c("gene1", "gene2", "gene3"))

  # Test with log transformation
  p_log <- getcountplot(
    df = mock_counts,
    intgroup = "condition",
    factor.levels = c("control", "treatment"),
    color = "gene",
    log = TRUE
  )
  expect_true(is.ggplot(p_log))

  # Test without log transformation
  p_linear <- getcountplot(
    df = mock_counts,
    intgroup = "condition",
    factor.levels = c("control", "treatment"),
    color = "gene",
    log = FALSE
  )
  expect_true(is.ggplot(p_linear))

  # Test with boxes disabled
  p_no_boxes <- getcountplot(
    df = mock_counts,
    intgroup = "condition",
    factor.levels = c("control", "treatment"),
    color = "gene",
    boxes = FALSE
  )
  expect_true(is.ggplot(p_no_boxes))

  # Test with legend disabled
  p_no_legend <- getcountplot(
    df = mock_counts,
    intgroup = "condition",
    factor.levels = c("control", "treatment"),
    color = "gene",
    legend = FALSE
  )
  expect_true(is.ggplot(p_no_legend))
})

test_that("getcountplot handles faceting correctly", {
  # Create mock data with batch information for faceting
  mock_counts <- create_mock_gene_counts(n_samples = 8, genes = c("gene1", "gene2"))

  # Test with single facet variable
  p_facet <- getcountplot(
    df = mock_counts,
    intgroup = "condition",
    factor.levels = c("control", "treatment"),
    color = "gene",
    facet = "batch"
  )
  expect_true(is.ggplot(p_facet))

  # Test with multiple facet variables
  p_multi_facet <- getcountplot(
    df = mock_counts,
    intgroup = "condition",
    factor.levels = c("control", "treatment"),
    color = "gene",
    facet = c("batch", "gene"),
    nrow = 2
  )
  expect_true(is.ggplot(p_multi_facet))

  # Test with free y-axis scales
  p_freey <- getcountplot(
    df = mock_counts,
    intgroup = "condition",
    factor.levels = c("control", "treatment"),
    color = "gene",
    facet = "batch",
    freey = TRUE
  )
  expect_true(is.ggplot(p_freey))
})

test_that("getcountplot handles trendlines correctly", {
  # Create mock data
  mock_counts <- create_mock_gene_counts(n_samples = 6, genes = c("gene1", "gene2"))

  # Test with smooth trendline (default)
  p_smooth <- getcountplot(
    df = mock_counts,
    intgroup = "condition",
    factor.levels = c("control", "treatment"),
    color = "gene",
    trendline = "smooth"
  )
  expect_true(is.ggplot(p_smooth))

  # Test with line trendline (median)
  p_line <- getcountplot(
    df = mock_counts,
    intgroup = "condition",
    factor.levels = c("control", "treatment"),
    color = "gene",
    trendline = "line"
  )
  expect_true(is.ggplot(p_line))

  # Test with no trendline
  p_no_trend <- getcountplot(
    df = mock_counts,
    intgroup = "condition",
    factor.levels = c("control", "treatment"),
    color = "gene",
    trendline = "none"
  )
  expect_true(is.ggplot(p_no_trend))
})

test_that("getcountplot handles axis limits and scaling correctly", {
  # Create mock data with known count ranges
  mock_counts <- create_mock_gene_counts(n_samples = 6, genes = c("gene1", "gene2"))

  # Test with custom y-axis limits
  p_limits <- getcountplot(
    df = mock_counts,
    intgroup = "condition",
    factor.levels = c("control", "treatment"),
    color = "gene",
    ymin = 50,
    ymax = 300,
    log = FALSE
  )
  expect_true(is.ggplot(p_limits))

  # Test with log scale and limits
  p_log_limits <- getcountplot(
    df = mock_counts,
    intgroup = "condition",
    factor.levels = c("control", "treatment"),
    color = "gene",
    ymin = 10,
    ymax = 1000,
    log = TRUE
  )
  expect_true(is.ggplot(p_log_limits))

  # Test with custom y-axis label
  p_custom_ylab <- getcountplot(
    df = mock_counts,
    intgroup = "condition",
    factor.levels = c("control", "treatment"),
    color = "gene",
    ylab = "Custom Y Label"
  )
  expect_true(is.ggplot(p_custom_ylab))
  expect_equal(p_custom_ylab$labels$y, "Custom Y Label")
})

test_that("getcountplot handles factor level filtering correctly", {
  # Create mock data with multiple conditions
  mock_counts <- data.frame(
    count = c(100, 120, 110, 200, 220, 210, 150, 160, 140),
    gene = rep(c("gene1", "gene2", "gene3"), each = 3),
    condition = rep(c("control", "treatment", "high_dose"), times = 3),
    sample = paste0("sample", 1:9),
    stringsAsFactors = FALSE
  )

  # Test filtering to specific factor levels
  p_filtered <- getcountplot(
    df = mock_counts,
    intgroup = "condition",
    factor.levels = c("control", "treatment"),  # Exclude "high_dose"
    color = "gene"
  )
  expect_true(is.ggplot(p_filtered))

  # Test with single factor level
  p_single <- getcountplot(
    df = mock_counts,
    intgroup = "condition",
    factor.levels = c("control"),
    color = "gene"
  )
  expect_true(is.ggplot(p_single))
})

test_that("getcountplot handles edge cases correctly", {
  # Test with minimal data
  minimal_counts <- data.frame(
    count = c(50, 60),
    gene = c("gene1", "gene1"),
    condition = c("control", "treatment"),
    sample = c("sample1", "sample2"),
    stringsAsFactors = FALSE
  )

  p_minimal <- getcountplot(
    df = minimal_counts,
    intgroup = "condition",
    factor.levels = c("control", "treatment"),
    color = "gene"
  )
  expect_true(is.ggplot(p_minimal))

  # Test with single gene
  single_gene <- create_mock_gene_counts(n_samples = 4, genes = "gene1")

  p_single_gene <- getcountplot(
    df = single_gene,
    intgroup = "condition",
    factor.levels = c("control", "treatment"),
    color = "gene"
  )
  expect_true(is.ggplot(p_single_gene))
})

# Test get_gene_counts function (used to generate data for getcountplot)
test_that("get_gene_counts extracts gene counts correctly", {
  # Create mock DESeq2 data
  mock_dds <- create_mock_dds(n_genes = 10, n_samples = 6)

  # Test basic gene count extraction
  gene_counts <- get_gene_counts(
    dds = mock_dds,
    gene = c("gene1", "gene2"),
    intgroup = "condition",
    norm_method = "libsize"
  )

  # Test that the result is a data frame
  expect_true(is.data.frame(gene_counts))

  # Test that it has the expected columns
  expected_cols <- c("count", "gene", "condition", "sample")
  expect_true(all(expected_cols %in% colnames(gene_counts)))

  # Test that it has the right number of rows (2 genes × 6 samples = 12 rows)
  expect_equal(nrow(gene_counts), 12)

  # Test that gene names are correct
  expect_true(all(c("gene1", "gene2") %in% gene_counts$gene))

  # Test that condition levels are correct
  expect_true(all(c("control", "treatment") %in% gene_counts$condition))
})

test_that("get_gene_counts handles different normalization methods", {
  # Create mock DESeq2 data
  mock_dds <- create_mock_dds(n_genes = 5, n_samples = 4)

  # Test with library size normalization
  counts_libsize <- get_gene_counts(
    dds = mock_dds,
    gene = "gene1",
    intgroup = "condition",
    norm_method = "libsize"
  )
  expect_true(is.data.frame(counts_libsize))
  expect_equal(nrow(counts_libsize), 4)  # 4 samples

  # Test with VST normalization
  counts_vst <- get_gene_counts(
    dds = mock_dds,
    gene = "gene1",
    intgroup = "condition",
    norm_method = "vst"
  )
  expect_true(is.data.frame(counts_vst))
  expect_equal(nrow(counts_vst), 4)  # 4 samples

  # The counts should be different between methods
  expect_false(identical(counts_libsize$count, counts_vst$count))
})

test_that("get_gene_counts handles missing genes gracefully", {
  # Create mock DESeq2 data
  mock_dds <- create_mock_dds(n_genes = 5, n_samples = 4)

  # Test with mix of existing and non-existing genes
  # Should skip non-existing genes and return data for existing ones
  gene_counts <- get_gene_counts(
    dds = mock_dds,
    gene = c("gene1", "nonexistent_gene", "gene2"),
    intgroup = "condition"
  )

  # Should only have data for gene1 and gene2
  expect_true(is.data.frame(gene_counts))
  expect_equal(nrow(gene_counts), 8)  # 2 genes × 4 samples = 8 rows
  expect_true(all(c("gene1", "gene2") %in% gene_counts$gene))
  expect_false("nonexistent_gene" %in% gene_counts$gene)
})

test_that("get_gene_counts validates intgroup parameter", {
  # Create mock DESeq2 data
  mock_dds <- create_mock_dds(n_genes = 3, n_samples = 4)

  # Test with valid intgroup
  expect_silent({
    gene_counts <- get_gene_counts(
      dds = mock_dds,
      gene = "gene1",
      intgroup = "condition"
    )
  })

  # Test with invalid intgroup
  expect_error(
    get_gene_counts(
      dds = mock_dds,
      gene = "gene1",
      intgroup = "nonexistent_column"
    ),
    "all variables in 'intgroup' must be columns of colData"
  )
})

# Integration test: get_gene_counts + getcountplot workflow
test_that("get_gene_counts and getcountplot work together correctly", {
  # Create mock DESeq2 data
  mock_dds <- create_mock_dds(n_genes = 5, n_samples = 6)

  # Extract gene counts using get_gene_counts
  gene_counts <- get_gene_counts(
    dds = mock_dds,
    gene = c("gene1", "gene2", "gene3"),
    intgroup = "condition",
    norm_method = "libsize"
  )

  # Use the extracted counts to create a plot with getcountplot
  p <- getcountplot(
    df = gene_counts,
    intgroup = "condition",
    factor.levels = c("control", "treatment"),
    color = "gene",
    ylab = "Normalized counts",
    title = "Gene Expression Plot"
  )

  # Test that the workflow produces a valid plot
  expect_true(is.ggplot(p))
  expect_equal(p$labels$title, "Gene Expression Plot")
  expect_equal(p$labels$y, "Normalized counts")

  # Test that the plot contains data for all 3 genes
  plot_data <- p$data
  expect_true(all(c("gene1", "gene2", "gene3") %in% plot_data$gene))

  # Test that the plot has the expected number of data points
  # 3 genes × 6 samples = 18 data points
  expect_equal(nrow(plot_data), 18)
})

# Test getcountplot with real DESeq2 workflow data
test_that("getcountplot works with DESeq2 transformed data", {
  # Create mock DESeq2 data and transform it
  mock_dds <- create_mock_dds()
  mock_rld <- varianceStabilizingTransformation(mock_dds, blind = TRUE)

  # Extract counts from transformed data
  vst_counts <- get_gene_counts(
    dds = mock_rld,
    gene = c("gene1", "gene2"),
    intgroup = "condition",
    norm_method = "vst"
  )

  # Create plot with VST-transformed data
  p_vst <- getcountplot(
    df = vst_counts,
    intgroup = "condition",
    factor.levels = c("control", "treatment"),
    color = "gene",
    log = FALSE,  # VST data is already log-like transformed
    ylab = "VST-transformed counts"
  )

  expect_true(is.ggplot(p_vst))
  expect_equal(p_vst$labels$y, "VST-transformed counts")
})

