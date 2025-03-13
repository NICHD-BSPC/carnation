#' Create object to be used with dashboard app
#'
#' @param config path to object creation config yaml
#'
#' @export
create_object <- function(config){
  config <- read_yaml(config)

  # - check if yaml blocks contain the allowed keys
  # - check if de_analysis 'counts' match counts 'key'
  # - check if paths exist
  check_config(config)

  # required keys
  counts <- config$counts
  de_analysis <- config$de_analysis

  enrich <- config$functional_enrichment
  degpatterns <- config$pattern_analysis

  # get annotation df
  annotation <- config$annotation
  id_col <- annotation$id
  symbol_col <- annotation$symbol
  id2symbol <- read.table(annotation$id2symbol,
                          sep='\t', header=TRUE)
  if(!all(c(id_col, symbol_col) %in% colnames(id2symbol))){
    stop(paste(id_col, 'and', symbol_col, 'columns not found in "id2symbol" file!'))
  }
  # replace NA symbols with id
  na_idx <- is.na(id2symbol[, symbol_col])
  id2symbol[na_idx, symbol_col] <- id2symbol[na_idx, id_col]
  rownames(id2symbol) <- id2symbol[, id_col]

  # make dds_list
  cat('Making list of DESeqDataSet\n')
  dds_list <- list()
  for(name in names(counts)){
    key <- counts[[name]][['key']]
    path <- counts[[name]][['path']]
    metadata <- counts[[name]][['metadata']]
    is_subset <- counts[[name]][['subset']]
    cat('\t', key, '\n')

    mat <- read.table(path, sep='\t',
                      header=TRUE, row.names=1)
    coldata <- read.table(metadata, sep='\t',
                          header=TRUE, row.names=1)

    # if subset flag is set, subset counts based on metadata
    if(is_subset){
      if(!all(rownames(coldata) %in% colnames(mat))){
        stop('All row names of metadata must be column names of counts matrix')
      }

      # keep counts from samples present in metadata
      # and implicitly reorder to match metadata row order
      mat <- mat[, rownames(coldata)]
    }

    # check that dimensions match
    if(!ncol(mat) == nrow(coldata)){
      stop('Count matrix must have the same number of columns as rows of the metadata file!')
    }

    # check that column names of count matrix and
    # rownames of metadata match
    if(!all(colnames(mat) == rownames(coldata))){
      stop('Column names of count matrix must match row names of metadata!')
    }

    # create DESeqDataSet object
    dds <- DESeqDataSetFromMatrix(mat,
                                  colData=coldata,
                                  design=~1)

    dds_list[[key]] <- dds
  }

  # make normalized rld_list
  cat('Normalizing count matrices\n')
  rld_list <- lapply(dds_list,
                  function(x)
                    varianceStabilizingTransformation(x, blind=TRUE))

  # make res_list
  cat('Making list of DE analysis results\n')
  res_list <- list()
  for(name in names(de_analysis)){
    key <- de_analysis[[name]][['key']]

    cat('\t', key, '\n')

    # read DE results and check column names from settings config
    res <- read.table(de_analysis[[name]][['path']], sep='\t',
                      header=TRUE, row.names=1)
    res <- check_res_columns(res, config, id2symbol)

    # check to make sure counts matrix and res object
    # have same genes
    counts_name <- de_analysis[[name]][['counts']]
    if(!all(rownames(dds_list[[counts_name]]) %in% rownames(res))){
      stop(paste0('Genes in count matrix: "', counts_name,
                  '" do not match genes in DE results: "', key,'"'))
    }

    res_list[[key]] <- list(
          'res'=res,
          'dds'=de_analysis[[name]][['counts']],
          'label'=de_analysis[[name]][['description']]
          )
  }

  # make initial object
  obj <- list(dds_list=dds_list,
              res_list=res_list,
              rld_list=rld_list)

  # make functional enrichment list
  if('functional_enrichment' %in% names(config)){
    cat('Making list of functional enrichment results\n')
    enrich_list <- list()
    genetonic <- list()
    for(name in names(enrich)){
      cat('\t', key, '\n')

      key <- enrich[[name]][['key']]
      res <- enrich[[name]][['res']]
      effect_class <- enrich[[name]][['effect_class']]
      ont <- enrich[[name]][['ontology']]

      # check to see if the 3 tiered list is created
      if(!is.list(enrich_list[[key]])){
        enrich_list[[key]] <- list()
        genetonic[[key]] <- list()
      }
      if(!is.list(enrich_list[[key]][[effect_class]])){
        enrich_list[[key]][[effect_class]] <- list()
        genetonic[[key]][[effect_class]] <- list()
      }

      # TODO: check to see if these are 'enrichResult' objects
      enrich_obj <- readRDS(enrich[[name]][['path']])
      enrich_list[[key]][[effect_class]][[ont]] <- enrich_obj

      # make genetonic object
      genetonic_obj <- enrich_to_genetonic(enrich_obj,
                              obj[['res_list']][[key]][['res']])
      genetonic[[key]][[effect_class]][[ont]] <- genetonic_obj
    }
    obj[['enrich_list']] <- enrich_list
    obj[['genetonic']] <- genetonic
  }

  # make pattern analysis list
  if('pattern_analysis' %in% names(config)){
    cat('Making list of pattern analysis results\n')
    degpatterns_list <- list()
    for(name in names(degpatterns)){
      key <- degpatterns[[name]][['key']]
      cat('\t', key, '\n')
      degpatterns_list[[key]] <- readRDS(degpatterns[[name]][['path']])
    }
    obj[['degpatterns']] <- degpatterns_list
  }

  return(obj)
}

#' check column names of DE analysis results
#'
#' - Required columns are c('padj', 'log2FoldChange', 'symbol')
#'
#' @param res data.frame with DE analysis results
#' @param config named list from reading object creation yaml
#' @param id2symbol data.frame with gene id -> symbol mapping
#'
#' @export
check_res_columns <- function(res, config, id2symbol){
  # if 'gene' column is missing, rownames of res are assumed
  # to be gene names
  if(!'gene' %in% colnames(res)){
    res <- res %>%
      mutate(gene=rownames(res)) %>%
      select(.data$gene, everything())
  }

  # if 'symbol' column is missing use annotation from config
  # and add it
  if(!'symbol' %in% colnames(res)){
    res$symbol <- id2symbol[rownames(res), config$annotation$symbol]
  }

  # check for padj & fold-change columns and rename to
  # 'padj' & 'log2FoldChange'
  padj_col <- config$settings$de_analysis$padj
  foldchange_col <- config$settings$de_analysis$foldchange
  if(!all(c(padj_col, foldchange_col) %in% colnames(res))){
    stop(paste(padj_col, 'and', foldchange_col, 'not found in DE results table!'))
  }
  res <- res %>%
    rename(padj = !!padj_col) %>%
    rename(log2FoldChange = !!foldchange_col)

  return(res)
}

#' validate object creation yaml file
#'
#' - check if yaml contains required & optional keys
#' - check if yaml blocks contain the allowed keys
#' - check if de_analysis 'counts' match counts 'key'
#' - check if paths exist
#'
#' @param config named list from reading object creation yaml
#' @param rules yaml file containing rules for object creation
#'
#' @export
check_config <- function(
                    config,
                    rules=system.file('extdata', 'create.yaml',
                                      package='dashboard_app')){
  rules <- read_yaml(rules)
  format <- rules$format
  optional <- rules$optional
  lower_level <- rules$lower_level

  if(!all(setdiff(names(format), optional) %in% names(config))){
    stop(paste0('Config yaml missing format keys: ',
                paste(setdiff(names(format), names(config)),
                      collapse=',')))
  }

  for(key in names(format)){
    for(key2 in names(format[[key]])){
      fun <- format[[key]][[key2]]

      if(!is.list(fun)){
        if(!key %in% lower_level){
          query <- config[[key]][[key2]]

          expr <- paste0(fun, '(query)')
          if(!eval(parse(text=expr))){
            stop(paste0('config:', key, ':', key2,
                        ' has incorrect format'))
          }
        } else {
          for(s in names(config[[key]])){
            query <- config[[key]][[s]][[key2]]

            expr <- paste0(fun, '(query)')
            if(!eval(parse(text=expr))){
              stop(paste0('config:', key, ':', s, ':', key2,
                          ' has incorrect format'))
            }
          }
        }
      } else {  # settings
        for(key3 in names(format[[key]][[key2]])){
          fun <- format[[key]][[key2]][[key3]]
          query <- config[[key]][[key2]][[key3]]


          expr <- paste0(fun, '(query)')
          if(!eval(parse(text=expr))){
            stop(paste0('config:', key, ':', key2, ':', 'key3',
                        ' has incorrect format'))
          }
        }

      }
    }
  }

  # 1. check de_analysis:counts values and pattern_analysis:counts
  #    values against counts:key
  # 2. check functional_analysis:res values against de_analysis:key

  # get all counts keys
  all_counts_keys <- NULL
  all_res_keys <- NULL
  for(name in names(config[['counts']])){
    all_counts_keys <- c(all_counts_keys,
                         config[['counts']][[name]][['key']])
  }

  for(name in names(config[['de_analysis']])){
    all_res_keys <- c(all_res_keys,
                     config[['de_analysis']][[name]][['key']])

    de_counts_key <- config[['de_analysis']][[name]][['counts']]
    if(!de_counts_key %in% all_counts_keys){
      stop(paste0('de_analysis counts key: "', de_counts_key,
                  '" not found in counts keys'))
    }
  }

  #for(name in names(config[['pattern_analysis']])){
  #  patterns_counts_key <- config[['pattern_analysis']][[name]][['counts']]
  #  if(!patterns_counts_key %in% all_counts_keys){
  #    stop(paste0('pattern_analysis counts key: "', patterns_counts_key,
  #                '" not found in counts keys'))
  #  }
  #}

  # check functional_analysis:res against de_analysis:keys
  for(name in names(config[['functional_analysis']])){
    functional_res_key <- config[['functional_analysis']][[name]][['res']]
    if(!functional_res_key %in% all_res_keys){
      stop(paste0('functional_analysis res key: "', functional_res_key,
                  '" not found in de_analysis keys'))
    }
  }

  cat('check_config: All checks passed\n')
}
