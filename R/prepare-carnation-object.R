#' Validate a carnation object
#'
#' This function takes various input data types (DE results, counts, enrichment,
#' pattern analysis) and validates them according to carnation's requirements,
#' returning a normalized intermediate object. Expensive derived-object
#' creation steps such as variance-stabilized counts and GeneTonic conversion
#' are handled separately by \code{materialize_carnation_object()}.
#'
#' @param res_list Named list of DE results. Each element should be either:
#'   \itemize{
#'     \item A data frame with DE results containing gene, symbol, pvalue, padj,
#'           log2FoldChange, and baseMean columns (or tool-specific alternatives)
#'     \item A list with slots: \code{res} (data frame), \code{dds} (name
#'           reference to dds_list element), \code{label} (comparison label)
#'   }
#' @param dds_list Named list of count data. Each element should be either:
#'   \itemize{
#'     \item A \code{DESeqDataSet} object
#'     \item A data frame or matrix of raw counts (first column=gene IDs, remaining=samples)
#'   }
#' @param rld_list Optional named list of variance-stabilized count objects.
#'   If NULL, these can be generated later via
#'   \code{materialize_carnation_object()}.
#' @param labels Optional named list of comparison labels. If NULL and res_list
#'   contains nested structure with \code{label} slots, labels will be extracted.
#' @param enrich_list Optional named list of functional enrichment results.
#'   Should be structured as: \code{enrich_list[[func_id]][[effect]][[pathway]]}.
#'   Each enrichment result must be a data frame in clusterProfiler format:
#'   \itemize{
#'     \item Over-representation: ID, Description, GeneRatio, BgRatio, pvalue,
#'           p.adjust, qvalue, geneID, Count
#'     \item GSEA: ID, Description, core_enrichment, setSize, pvalue, p.adjust,
#'           qvalue, NES
#'   }
#' @param degpatterns Optional named list of pattern analysis results. Each
#'   element should be either a data frame or a list with \code{$normalized}
#'   slot containing a data frame with columns: genes, value, and either cluster
#'   or columns starting with "cutoff".
#' @param metadata Optional data frame with sample metadata. Required if
#'   \code{dds_list} contains count matrices instead of DESeqDataSet objects.
#'   First column should be sample names matching column names in count matrices.
#' @param dds_mapping Optional named list mapping \code{res_list} elements to
#'   \code{dds_list} objects. Required if \code{res_list} is a list of data frames.
#' @param config Optional config list. If NULL, will use \code{get_config()},
#'   including any supported local config overrides.
#'
#' @return A validated list with canonical slots \code{res_list},
#'   \code{dds_list}, optional \code{rld_list}, \code{labels},
#'   \code{dds_mapping}, \code{enrich_list}, \code{degpatterns}, and
#'   \code{metadata} when supplied.
#'
#' @details
#' This function performs comprehensive validation of all input data:
#' \itemize{
#'   \item DE results: Checks for required columns (with support for DESeq2,
#'         edgeR, limma), ensures gene and symbol columns exist
#'   \item Counts: Validates structure, checks sample name matching with metadata
#'   \item Enrichment: Validates clusterProfiler format (OR or GSEA)
#'   \item Pattern analysis: Checks for required columns (genes, value, cluster)
#' }
#'
#' If validation fails, the function will stop with an informative error message.
#'
#' @examples
#' \dontrun{
#' # Minimal example with DE results and counts
#' library(DESeq2)
#'
#' # Create example data
#' dds <- makeExampleDESeqDataSet()
#' dds <- DESeq(dds)
#' res <- results(dds, contrast = c("condition", "A", "B"))
#' rld <- varianceStabilizingTransformation(dds, blind = TRUE)
#'
#' # Validate object inputs
#' obj <- validate_carnation_object(
#'   res_list = list(
#'     comp1 = list(
#'       res = as.data.frame(res),
#'       dds = "main",
#'       label = "A vs B"
#'     )
#'   ),
#'   dds_list = list(main = dds),
#'   rld_list = list(main = rld)
#' )
#'
#' materialized <- materialize_carnation_object(obj, cores = 1)
#' final_obj <- make_final_object(materialized)
#'
#' # Save for use with carnation
#' saveRDS(final_obj, "my_analysis.rds")
#'
#' # Alternative: start from count matrix and metadata
#' counts <- as.data.frame(counts(dds))
#' counts$gene <- rownames(counts)
#' counts <- counts[, c(ncol(counts), 1:(ncol(counts)-1))]
#' metadata <- as.data.frame(colData(dds))
#' metadata$sample <- rownames(metadata)
#' metadata <- metadata[, c(ncol(metadata), 1:(ncol(metadata)-1))]
#'
#' obj <- validate_carnation_object(
#'   res_list = list(comp1 = as.data.frame(res)),
#'   dds_list = list(main = counts),
#'   metadata = metadata,
#'   dds_mapping = list(comp1 = "main")
#' )
#' }
#'
#' @rdname validate_carnation_object
#' @name validate_carnation_object
#'
#' @return A list containing normalized inputs with elements
#'   \code{res_list}, \code{dds_list}, optional \code{rld_list},
#'   \code{labels}, \code{dds_mapping}, and optional \code{enrich_list},
#'   \code{degpatterns}, and \code{metadata}.
#'
#' @export
validate_carnation_object <- function(res_list,
                                      dds_list,
                                      rld_list = NULL,
                                      labels = NULL,
                                      enrich_list = NULL,
                                      degpatterns = NULL,
                                      metadata = NULL,
                                      dds_mapping = NULL,
                                      config = NULL) {
  if (is.null(config)) {
    config <- get_config()
  }

  # helper to check if arg is a named list without duplicated names
  .validate_named_input <- function(x, arg) {
    if (!is.list(x) || is.null(names(x))) {
      stop(arg, " must be a named list")
    }

    nms <- names(x)
    if (any(is.na(nms)) || any(nms == "")) {
      stop(arg, " must have non-empty names")
    }
    if (any(duplicated(nms))) {
      dupes <- unique(nms[duplicated(nms)])
      stop(arg, " contains duplicate names: ",
           paste(dupes, collapse = ", "))
    }
    if (any(grepl("[[:space:],]", nms))) {
      bad <- unique(nms[grepl("[[:space:],]", nms)])
      stop(arg, " names cannot contain white-space or commas: ",
           paste(bad, collapse = ", "))
    }
  }

  obj <- list()

  if (missing(res_list) || is.null(res_list)) {
    stop("res_list is required")
  }
  .validate_named_input(res_list, "res_list")

  if (missing(dds_list) || is.null(dds_list)) {
    stop("dds_list is required")
  }
  .validate_named_input(dds_list, "dds_list")

  if (!is.null(rld_list)) {
    .validate_named_input(rld_list, "rld_list")
  }

  # if 'labels' is specified, make sure there are no missing ones
  if (!is.null(labels)) {
    .validate_named_input(labels, "labels")
    if (!all(names(res_list) %in% names(labels))) {
      missing_labels <- setdiff(names(res_list), names(labels))
      stop("labels is missing entries for: ",
           paste(missing_labels, collapse = ", "))
    }
  }

  # if dds_mapping is specified, check to see if all dds_list/res_list
  # elements referenced actually exist
  if(!is.null(dds_mapping)){
    .validate_named_input(dds_mapping, "dds_mapping")
    res_names <- names(dds_mapping)
    dds_names <- unname(unique(unlist(dds_mapping)))
    if(!all(res_names %in% names(res_list))){
      stop("Some dds_mapping keys do not exist in res_list:",
           paste(setdiff(res_names, names(res_list)), collapse=','))
    } else if(!(all(dds_names %in% names(dds_list)))){
      stop("Some dds_mapping values do not exist in dds_list:",
           paste(setdiff(dds_names, names(dds_list)), collapse=','))
    }
  }

  message("Validating DE results...")

  validated_res_list <- list()
  for (res_name in names(res_list)) {
    res_item <- res_list[[res_name]]

    if (is.list(res_item) && all(c("res", "dds", "label") %in% names(res_item))) {
      res <- res_item$res
      dds_ref <- res_item$dds
      label <- if (!is.null(labels)) labels[[res_name]] else res_item$label
    } else if (is.data.frame(res_item)) {
      res <- res_item
      if(length(dds_list) == 1){
        dds_ref <- names(dds_list)[1]
      } else if(!is.null(dds_mapping)){
        dds_ref <- dds_mapping[[ res_name ]]
      } else {
        stop("If res_list is list of data frames and length(dds_list) > 1,",
             " dds_mapping must be specified")
      }
      label <- if (!is.null(labels)) labels[[res_name]] else res_name
    } else {
      stop("Each element in res_list must be a data frame or a list with slots: res, dds, label")
    }

    res <- as.data.frame(res)

    if (!is.character(dds_ref) || length(dds_ref) != 1 || is.na(dds_ref) ||
        !nzchar(dds_ref)) {
      stop("DE result '", res_name,
           "' must reference a single dds_list name in its `dds` field")
    }
    if (!dds_ref %in% names(dds_list)) {
      stop("DE result '", res_name, "' references unknown dds_list entry: ",
           dds_ref)
    }
    if (!is.character(label) || length(label) != 1 || is.na(label) ||
        !nzchar(label)) {
      stop("DE result '", res_name, "' must have a non-empty label")
    }

    gene_idx <- which(toupper(colnames(res)) %in% "GENE")
    if (length(gene_idx) == 0 && is.null(rownames(res))) {
      stop("DE results '", res_name, "' must have a 'gene' column or row names")
    }

    if (length(gene_idx) == 0) {
      res$gene <- rownames(res)
    }

    # check for symbols
    # 1. plug in gene names if symbol column not found
    # 2. handle NAs by using gene names/rownames instead
    symbol_idx <- which(toupper(colnames(res)) %in% "SYMBOL")
    gene_idx <- which(toupper(colnames(res)) %in% "GENE")
    if (length(symbol_idx) == 0) {
      message("  Adding 'symbol' column from 'gene' column for: ", res_name)
      if (length(gene_idx) > 0) {
        res$symbol <- res[, gene_idx[1]]
      } else {
        res$symbol <- rownames(res)
      }
    } else if (any(is.na(res[, symbol_idx[1]]))) {
      na_symbol_idx <- is.na(res[, symbol_idx[1]])
      if (length(gene_idx) > 0) {
        res[na_symbol_idx, symbol_idx[1]] <-
          res[na_symbol_idx, gene_idx[1]]
      } else if (!is.null(rownames(res))) {
        res[na_symbol_idx, symbol_idx[1]] <-
          rownames(res)[na_symbol_idx]
      }
    }

    # check for supported columns
    res_validated <- .check_res_columns(res, config$server$de_analysis$column_names)
    if (length(res_validated$missing_cols) > 0) {
      stop("DE results '", res_name, "' is missing required columns: ",
           paste(res_validated$missing_cols, collapse = ", "),
           "\nSupported alternatives:\n",
           paste(capture.output(str(config$server$de_analysis$column_names)),
                 collapse = "\n"))
    }

    # order by padj then plug back in
    validated_res_list[[res_name]] <- list(
      res = res_validated$res[order(res_validated$res$padj), ],
      dds = dds_ref,
      label = label
    )

    message("  ✓ ", res_name)
  }

  obj$res_list <- validated_res_list

  message("\nValidating count data...")

  validated_dds_list <- list()
  validated_rld_list <- list()
  metadata_df <- NULL

  if (!is.null(metadata)) {
    metadata_df <- as.data.frame(metadata)
    if (ncol(metadata_df) < 1) {
      stop("metadata must have at least one column (sample names)")
    }
    if (any(duplicated(metadata_df[, 1]))) {
      stop("metadata contains duplicate sample names in first column")
    }
    rownames(metadata_df) <- metadata_df[, 1]
  }

  for (dds_name in names(dds_list)) {
    dds_item <- dds_list[[dds_name]]
    rld_item <- NULL
    if (!is.null(rld_list) && dds_name %in% names(rld_list)) {
      rld_item <- rld_list[[dds_name]]
    }

    if (inherits(dds_item, "DESeqDataSet")) {
      validated_dds_list[[dds_name]] <- dds_item

      if (!is.null(rld_item)) {
        if (!inherits(rld_item, "DESeqTransform")) {
          stop("rld_list entry '", dds_name, "' must be a DESeqTransform object")
        }
        if (!identical(rownames(rld_item), rownames(dds_item))) {
          stop("rld_list entry '", dds_name,
               "' must have the same row names as its dds_list entry")
        }
        if (!identical(colnames(rld_item), colnames(dds_item))) {
          stop("rld_list entry '", dds_name,
               "' must have the same sample columns as its dds_list entry")
        }
        validated_rld_list[[dds_name]] <- rld_item
      }
    } else if (is.data.frame(dds_item) || is.matrix(dds_item)) {
      if (is.null(metadata_df)) {
        stop("metadata is required when dds_list contains count matrices (for: ",
             dds_name, ")")
      }

      counts_df <- as.data.frame(dds_item)
      if (ncol(counts_df) < 2) {
        stop("Count matrix '", dds_name,
             "' must contain a gene ID column and at least one sample column")
      }
      if (!is.character(counts_df[, 1])) {
        stop("First column of count matrix '", dds_name,
             "' must be character (gene IDs)")
      }
      if (any(duplicated(counts_df[, 1]))) {
        stop("Count matrix '", dds_name, "' contains duplicate gene IDs")
      }

      count_cols <- 2:ncol(counts_df)
      if (!all(apply(counts_df[, count_cols], 2, is.numeric))) {
        stop("Columns 2+ in count matrix '", dds_name, "' must all be numeric")
      }

      gene_ids <- counts_df[, 1]
      counts_only <- counts_df[, -1, drop = FALSE]

      if (nrow(metadata_df) != ncol(counts_only)) {
        stop("Number of samples in metadata (", nrow(metadata_df),
             ") does not match number of columns in count matrix '", dds_name,
             "' (", ncol(counts_only), ")")
      }
      if (!all(colnames(counts_only) %in% rownames(metadata_df))) {
        missing <- setdiff(colnames(counts_only), rownames(metadata_df))
        stop("Sample names in count matrix '", dds_name,
             "' not found in metadata: ", paste(missing, collapse = ", "))
      }
      if (!all(rownames(metadata_df) %in% colnames(counts_only))) {
        missing <- setdiff(rownames(metadata_df), colnames(counts_only))
        stop("Sample names in metadata not found in count matrix '", dds_name,
             "': ", paste(missing, collapse = ", "))
      }

      validated_dds_list[[dds_name]] <- counts_df

      if (!is.null(rld_item)) {
        if (!inherits(rld_item, "DESeqTransform")) {
          stop("rld_list entry '", dds_name, "' must be a DESeqTransform object")
        }
        if (!identical(rownames(rld_item), gene_ids)) {
          stop("rld_list entry '", dds_name,
               "' must have the same row names as its dds_list entry")
        }
        if (!identical(colnames(rld_item), colnames(counts_only))) {
          stop("rld_list entry '", dds_name,
               "' must have the same sample columns as its dds_list entry")
        }
        validated_rld_list[[dds_name]] <- rld_item
      }
    } else {
      stop("dds_list element '", dds_name,
           "' must be a DESeqDataSet or count matrix/data frame")
    }

    if (inherits(validated_dds_list[[dds_name]], "DESeqDataSet") &&
        is.null(SummarizedExperiment::colData(validated_dds_list[[dds_name]])$sample)) {
      SummarizedExperiment::colData(validated_dds_list[[dds_name]])$sample <-
        rownames(SummarizedExperiment::colData(validated_dds_list[[dds_name]]))
    }
    if (!is.null(validated_rld_list[[dds_name]]) &&
        is.null(SummarizedExperiment::colData(validated_rld_list[[dds_name]])$sample)) {
      SummarizedExperiment::colData(validated_rld_list[[dds_name]])$sample <-
        rownames(SummarizedExperiment::colData(validated_rld_list[[dds_name]]))
    }

    message("  ✓ ", dds_name)
  }

  obj$dds_list <- validated_dds_list
  if (length(validated_rld_list) > 0) {
    obj$rld_list <- validated_rld_list
  }
  if (!is.null(metadata_df)) {
    obj$metadata <- metadata_df
  }
  obj$labels <- if (is.null(labels)) {
    lapply(validated_res_list, function(x) x$label)
  } else {
    labels[names(validated_res_list)]
  }
  obj$dds_mapping <- lapply(validated_res_list, function(x) x$dds)

  if (!is.null(enrich_list)) {
    message("\nValidating functional enrichment results...")

    .validate_named_input(enrich_list, "enrich_list")

    validated_enrich_list <- list()

    cprof_or_cols <- c("ID", "Description", "GeneRatio", "BgRatio",
                       "pvalue", "p.adjust", "qvalue", "geneID", "Count")
    cprof_gsea_cols <- c("ID", "Description", "core_enrichment", "setSize",
                         "pvalue", "p.adjust", "qvalue", "NES")

    for (func_id in names(enrich_list)) {
      func_item <- enrich_list[[func_id]]

      if (!is.list(func_item)) {
        stop("enrich_list element '", func_id, "' must be a list")
      }

      validated_enrich_list[[func_id]] <- list()

      # check enrich_list elements
      # if 'res' key present, save or try to guess from
      # res_list names
      func_item_names <- names(func_item)
      res_ref <- NULL
      if ("res" %in% names(func_item)) {
        res_ref <- func_item$res
        func_item_names <- setdiff(func_item_names, 'res')
      } else if (func_id %in% names(validated_res_list)) {
        res_ref <- func_id
      } else if (length(validated_res_list) == 1) {
        res_ref <- names(validated_res_list)[1]
      }

      for (effect in func_item_names) {
        effect_item <- func_item[[effect]]
        if (!is.list(effect_item) || is.null(names(effect_item))) {
          stop("enrich_list[[", func_id, "]][[", effect, "]] must be a list")
        }
        if (any(is.na(names(effect_item))) || any(names(effect_item) == "")) {
          stop("enrich_list[[", func_id, "]][[", effect,
               "]] must have non-empty pathway names")
        }
        if (any(duplicated(names(effect_item)))) {
          stop("enrich_list[[", func_id, "]][[", effect,
               "]] contains duplicate pathway names")
        }

        validated_enrich_list[[func_id]][[effect]] <- list()

        for (pathway in names(effect_item)) {
          enrich_df <- effect_item[[pathway]]

          if (inherits(enrich_df, "enrichResult") || inherits(enrich_df, "gseaResult")) {
            enrich_df <- enrich_df@result
          } else if (!is.data.frame(enrich_df)) {
            stop("enrich_list[[", func_id, "]][[", effect, "]][[", pathway,
                 "]] must be a data frame or enrichResult/gseaResult object")
          }

          if (all(cprof_or_cols %in% colnames(enrich_df))) {
            rownames(enrich_df) <- enrich_df$ID
          } else if (all(cprof_gsea_cols %in% colnames(enrich_df))) {
            rownames(enrich_df) <- enrich_df$ID
          } else {
            stop("Enrichment result '", func_id, ".", effect, ".", pathway,
                 "' does not match clusterProfiler format.\n",
                 "Required columns for OR: ", paste(cprof_or_cols, collapse = ", "), "\n",
                 "Required columns for GSEA: ", paste(cprof_gsea_cols, collapse = ", "))
          }

          validated_enrich_list[[func_id]][[effect]][[pathway]] <- enrich_df
          if (!is.null(res_ref)) {
            if (!res_ref %in% names(validated_res_list)) {
              stop("enrich_list element '", func_id,
                   "' references unknown DE result: ", res_ref)
            }
            validated_enrich_list[[func_id]]$res <- res_ref
          }

          message("  ✓ ", func_id, ".", effect, ".", pathway)
        }
      }
    }

    obj$enrich_list <- validated_enrich_list
  }

  if (!is.null(degpatterns)) {
    message("\nValidating pattern analysis results...")

    .validate_named_input(degpatterns, "degpatterns")

    validated_degpatterns <- list()
    for (dp_name in names(degpatterns)) {
      dp_item <- degpatterns[[dp_name]]
      if (!is_valid_pattern_obj(dp_item, require_symbol = FALSE)) {
        stop("Pattern analysis object '", dp_name, "' failed validation")
      }

      if (is.list(dp_item) && "normalized" %in% names(dp_item)) {
        validated_degpatterns[[dp_name]] <- dp_item$normalized
      } else {
        validated_degpatterns[[dp_name]] <- dp_item
      }
      message("  ✓ ", dp_name)
    }

    obj$degpatterns <- validated_degpatterns
  }

  obj
}

#' Materialize expensive carnation object components
#'
#' This function materializes expensive derived pieces for a validated
#' carnation object, including DESeqDataSet creation from raw count matrices,
#' variance-stabilized counts, and GeneTonic conversions.
#'
#' @param obj A validated object returned by
#'   \code{validate_carnation_object()} or
#'   \code{validate_loaded_carnation_object()}.
#' @param config Optional config list. If NULL, will use \code{get_config()}.
#' @param cores Optional number of worker processes. If NULL, uses
#'   \code{config$server$cores}.
#'
#' @examplesIf interactive()
#' # Minimal example with DE results and counts
#' library(DESeq2)
#'
#' # Create example data
#' dds <- makeExampleDESeqDataSet()
#' dds <- DESeq(dds)
#' res <- results(dds, contrast = c("condition", "A", "B"))
#' rld <- varianceStabilizingTransformation(dds, blind = TRUE)
#'
#' # Validate object inputs
#' obj <- validate_carnation_object(
#'   res_list = list(
#'     comp1 = list(
#'       res = as.data.frame(res),
#'       dds = "main",
#'       label = "A vs B"
#'     )
#'   ),
#'   dds_list = list(main = dds),
#'   rld_list = list(main = rld)
#' )
#'
#' materialized <- materialize_carnation_object(obj, cores = 1)
#'
#' @return The input object with materialized \code{dds_list},
#'   \code{rld_list}, and optional \code{genetonic} slots.
#'
#' @export
materialize_carnation_object <- function(obj,
                                         config = NULL,
                                         cores = NULL) {
  if (is.null(config)) {
    config <- get_config()
  }
  if (is.null(cores)) {
    cores <- config$server$cores
  }

  if (!is.list(obj) || is.null(obj$res_list) || is.null(obj$dds_list)) {
    stop("obj must be a validated carnation object with res_list and dds_list")
  }

  metadata_df <- obj$metadata
  rld_list <- obj$rld_list
  if (is.null(rld_list)) {
    rld_list <- list()
  }

  dds_names <- names(obj$dds_list)
  materialized_counts <- BiocParallel::bplapply(
    dds_names,
    function(dds_name) {
      dds_item <- obj$dds_list[[dds_name]]
      rld_item <- rld_list[[dds_name]]

      if (inherits(dds_item, "DESeqDataSet")) {
        dds_obj <- dds_item
      } else if (is.data.frame(dds_item) || is.matrix(dds_item)) {
        if (is.null(metadata_df)) {
          stop("metadata is required when dds_list contains count matrices")
        }

        counts_df <- as.data.frame(dds_item)
        rownames(counts_df) <- counts_df[, 1]
        counts_df <- counts_df[, -1, drop = FALSE]
        counts_df <- counts_df[, rownames(metadata_df), drop = FALSE]

        dds_obj <- DESeq2::DESeqDataSetFromMatrix(
          countData = as.matrix(counts_df),
          colData = metadata_df,
          design = ~ 1
        )
      } else {
        stop("dds_list element '", dds_name,
             "' must be a DESeqDataSet or count matrix/data frame")
      }

      if (is.null(SummarizedExperiment::colData(dds_obj)$sample)) {
        SummarizedExperiment::colData(dds_obj)$sample <-
          rownames(SummarizedExperiment::colData(dds_obj))
      }

      if (is.null(rld_item)) {
        rld_obj <- DESeq2::varianceStabilizingTransformation(dds_obj, blind = TRUE)
      } else {
        rld_obj <- rld_item
      }

      if (!inherits(rld_obj, "DESeqTransform")) {
        stop("rld_list entry '", dds_name, "' must be a DESeqTransform object")
      }
      if (!identical(rownames(rld_obj), rownames(dds_obj))) {
        stop("rld_list entry '", dds_name,
             "' must have the same row names as its dds_list entry")
      }
      if (!identical(colnames(rld_obj), colnames(dds_obj))) {
        stop("rld_list entry '", dds_name,
             "' must have the same sample columns as its dds_list entry")
      }
      if (is.null(SummarizedExperiment::colData(rld_obj)$sample)) {
        SummarizedExperiment::colData(rld_obj)$sample <-
          rownames(SummarizedExperiment::colData(rld_obj))
      }

      list(dds = dds_obj, rld = rld_obj)
    },
    BPPARAM = BiocParallel::MulticoreParam(max(1, cores))
  )
  names(materialized_counts) <- dds_names

  obj$dds_list <- setNames(lapply(materialized_counts, function(x) x$dds), dds_names)
  obj$rld_list <- setNames(lapply(materialized_counts, function(x) x$rld), dds_names)

  if (!is.null(obj$enrich_list) && is.null(obj$genetonic)) {
    obj$genetonic <- .materialize_genetonic_list(
      enrich_list = obj$enrich_list,
      res_list = obj$res_list,
      cores = cores
    )
  }

  obj
}

# Internal function to validate loaded carnation object
validate_loaded_carnation_object <- function(obj,
                                             config = NULL) {
  if (is.null(config)) {
    config <- get_config()
  }

  if (!is.list(obj)) {
    stop("Loaded object must be a list")
  }

  n <- names(obj)
  res_name <- n[grep("res", n)]
  dds_name <- n[setdiff(grep("dds", n),
                        c(grep("all_dds", n),
                          grep("dds_mapping", n)))]
  rld_name <- n[setdiff(grep("rld", n), grep("all_rld", n))]
  enrich_name <- n[grep("enrich", n)]
  degpatterns_name <- n[grep("degpatterns", n)]

  if (length(res_name) == 0 || length(dds_name) == 0) {
    missing <- if (length(res_name) == 0 && length(dds_name) == 0) {
      "both"
    } else if (length(res_name) == 0) {
      "DE result"
    } else {
      "counts table"
    }
    stop("Loaded object must have both DE results & counts table. Missing: ",
         missing)
  }
  if (length(dds_name) != 1) {
    stop("Uploaded object must have exactly 1 counts list")
  }
  if (length(rld_name) > 1) {
    stop("Uploaded object must have exactly 1 normalized counts list")
  }

  labels <- obj$labels
  dds_mapping <- obj$dds_mapping

  validate_carnation_object(
    res_list = obj[[res_name]],
    dds_list = obj[[dds_name]],
    rld_list = if (length(rld_name) == 1) obj[[rld_name]] else NULL,
    labels = labels,
    enrich_list = if (length(enrich_name) > 0) obj[[enrich_name]] else NULL,
    degpatterns = if (length(degpatterns_name) > 0) obj[[degpatterns_name]] else NULL,
    metadata = obj$metadata,
    dds_mapping = dds_mapping,
    config = config
  )
}

.materialize_genetonic_list <- function(enrich_list,
                                        res_list,
                                        cores) {
  elem_names <- character()
  sep <- "*"
  res_keys <- list()

  for (x in names(enrich_list)) {
    for (y in names(enrich_list[[x]])) {
      if (y == "res") {
        res_keys[[x]] <- enrich_list[[x]][["res"]]
        next
      }
      for (z in names(enrich_list[[x]][[y]])) {
        elem_names <- c(elem_names, paste(x, y, z, sep = sep))
      }
    }
  }
  names(elem_names) <- elem_names

  flat_obj <- BiocParallel::bplapply(
    elem_names,
    function(x) {
      toks <- strsplit(x, split = sep, fixed = TRUE)[[1]]
      res_ref <- toks[1]
      if (!(res_ref %in% names(res_list)) && toks[1] %in% names(res_keys)) {
        res_ref <- res_keys[[toks[1]]]
      }

      eres <- enrich_list[[toks[1]]][[toks[2]]][[toks[3]]]
      if (all(c("ID", "Description", "core_enrichment", "setSize",
                "pvalue", "p.adjust", "qvalue", "NES") %in% colnames(eres))) {
        enrich_obj <- makeEnrichResult(eres, type = "gseaResult")
      } else {
        enrich_obj <- makeEnrichResult(eres, type = "enrichResult")
      }

      if (res_ref %in% names(res_list)) {
        enrich_to_genetonic(enrich_obj, res_list[[res_ref]]$res)
      } else {
        dummy_genetonic(enrich_obj)
      }
    },
    BPPARAM = BiocParallel::MulticoreParam(max(1, cores))
  )
  names(flat_obj) <- elem_names

  genetonic <- list()
  for (x in names(enrich_list)) {
    genetonic[[x]] <- list()
    for (y in names(enrich_list[[x]])) {
      if (y == "res") {
        next
      }
      genetonic[[x]][[y]] <- list()
      for (z in names(enrich_list[[x]][[y]])) {
        genetonic[[x]][[y]][[z]] <- flat_obj[[paste(x, y, z, sep = sep)]]
      }
    }
  }

  genetonic
}
