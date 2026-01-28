#' Create access yaml
#'
#' This function creates an access yaml file.
#' This is primarily intended for the first run.
#'
#' @param user User name
#' @param user_group User group
#' @param data_area Path to data area containing RDS files
#'
#' @return Invisibly returns `NULL`. This function is primarily used for
#'   its side effect of saving a yaml file with access settings
#'
#' @examples
#' # save access details to file
#' home <- Sys.getenv('HOME')
#'
#' # create carnation data area if it doesn't exist
#' carnation_home <- file.path(home, 'carnation/data')
#' if(!dir.exists(carnation_home)) dir.create(carnation_home)
#'
#' create_access_yaml(user = 'admin',
#'                    user_group = 'admin',
#'                    data_area = carnation_home)
#'
#' @export
create_access_yaml <- function(user, user_group, data_area){
  ug <- setNames(as.list(user_group), user)
  da <- setNames(as.list(data_area), user_group)

  path <- get_access_path()

  write_yaml(list(user_group=ug, data_area=da),
             path)
}

#' Read access yaml with user groups and data areas
#'
#' This function reads the access yaml file and
#' returns user groups and data areas
#' as a list of data frames.
#'
#' @return return carnation access settings from yaml file
#'
#' @examples
#' # save access details to file
#' home <- Sys.getenv('HOME')
#'
#' # create carnation data area if it doesn't exist
#' carnation_home <- file.path(home, 'carnation/data')
#' if(!dir.exists(carnation_home)) dir.create(carnation_home)
#'
#' create_access_yaml(user = 'admin',
#'                    user_group = 'admin',
#'                    data_area = carnation_home)
#'
#' al <- read_access_yaml()
#'
#' @export
read_access_yaml <- function(){
  # get path to access yaml
  f <- get_access_path()

  # check if yaml exists
  if(!file.exists(f)){
    stop('Access yaml not found. Have you run "create_access_yaml()" yet?')
  }

  al <- read_yaml(f)

  return(al)
}

#' Save access yaml to file
#'
#' This function saves access details (user groups
#' and data areas) to the designated access yaml file.
#'
#' @param lst list of data frames with user_groups and
#'  data_areas
#'
#' @return save access settings to yaml file
#'
#' @examples
#' # save access details to file
#' home <- Sys.getenv('HOME')
#'
#' # create carnation data area if it doesn't exist
#' carnation_home <- file.path(home, 'carnation/data')
#' if(!dir.exists(carnation_home)) dir.create(carnation_home)
#'
#' create_access_yaml(user = 'admin',
#'                    user_group = 'admin',
#'                    data_area = carnation_home)
#'
#' # read access yaml
#' lst <- read_access_yaml()
#'
#' # add new user
#' lst$user_group$admin <- c(lst$user_group$admin, 'user1')
#'
#' # save to access settings
#' save_access_yaml(lst)
#'
#' @export
save_access_yaml <- function(lst){
  # get access file
  f <- get_access_path()

  write_yaml(list(user_group=lst$user_group,
                  data_area=lst$data_area), f)
}

#' is user an admin?
#'
#' @param u username
#'
#' @return boolean to indicate is user is in admin group
#'
#' @examples
#' # check if default user is admin
#' yy <- is_site_admin(u='admin')
#'
#' @export
is_site_admin <- function(u){
  cfg <- get_config()
  admin <- cfg$server$site_admin

  if(is.null(u)) return(FALSE)
  if(u %in% admin) return(TRUE)
  else return(FALSE)
}

#' is user is in admin group?
#'
#' @param u username
#'
#' @return TRUE/FALSE to indicate if the user is part of the admin group
#'
#' @examples
#' # save access details to file
#' home <- Sys.getenv('HOME')
#'
#' # create carnation data area if it doesn't exist
#' carnation_home <- file.path(home, 'carnation/data')
#' if(!dir.exists(carnation_home)) dir.create(carnation_home)
#'
#' create_access_yaml(user = 'admin',
#'                    user_group = 'admin',
#'                    data_area = carnation_home)
#'
#' check <- in_admin_group('user')
#'
#' @export
in_admin_group <- function(u){
  al <- read_access_yaml()

  cfg <- get_config()
  admin_group <- cfg$server$admin_group

  ll <- check_user_access(al, u, admin_group)

  if(admin_group %in% ll$user_group) return(TRUE)
  else return(FALSE)
}

#' Get config
#'
#' This function reads the config.yaml and
#' returns the list
#'
#' @return list containing config items
#'
#' @examples
#' cfg <- get_config()
#'
#' @export
get_config <- function(){
  cfg_path <- system.file('extdata', 'config.yaml',
                          package=packageName())
  cfg <- read_yaml(cfg_path)
  cfg
}

#' Get project name from path
#'
#' This function takes in a path to an RDS file and returns
#' a string to be used as project name
#'
#' @param x character path to RDS file
#' @param depth integer how many levels below path to look?
#' @param end_offset integer how far from the end of path to end?
#' @param staging_dir name of staging directory
#' @param fsep file separator to split path with
#'
#' @return project name parsed from path to object
#'
#' @examples
#' # path to carnation object
#' obj_path <- "/path/to/project/test/main.rnaseq.rds"
#'
#' # parsed project name
#' get_project_name_from_path(obj_path, depth = 2, end_offset = 0)
#'
#' @export
get_project_name_from_path <- function(x,
                                       depth=2, end_offset=0,
                                       staging_dir='dev',
                                       fsep=.Platform$file.sep){
  d <- dirname(x)

  # split path by fsep
  tok <- unlist(strsplit(d, fsep))

  # if path contains staging_dir, increase depth
  if(any(tok == staging_dir)){
    depth <- depth + 1
  }

  # join upto 'level' elements from the end
  end_idx <- length(tok) - end_offset
  start_idx <- length(tok) - depth + 1

  p <- paste(tok[start_idx:end_idx], collapse=fsep)

  p
}

#' Get read counts for gene
#'
#' This is a simple function to obtain read counts for a
#' specified gene, based on the DESeq2::plotCounts function.
#'
#' @param dds DESeqDataSet object
#' @param gene gene name vector
#' @param intgroup metadata variable to attach to counts
#' @param norm_method normalization method, can be 'libsize' (default) or 'vst'
#'
#' @return data.frame with gene counts
#'
#' @examples
#' # make example DESeq data set
#' dds <- DESeq2::makeExampleDESeqDataSet()
#'
#' # get counts for gene1
#' gg <- get_gene_counts(dds, 'gene1')
#'
#' @export
get_gene_counts <- function (dds,
                             gene,
                             intgroup = "condition",
                             norm_method = 'libsize')
{

  if(norm_method == 'libsize'){
    normalized <-TRUE
  } else if(norm_method == 'vst'){
    normalized <- FALSE
  }

  if (!all(intgroup %in% names(colData(dds))))
    stop("all variables in 'intgroup' must be columns of colData")

  if(inherits(dds, 'DESeqDataSet')){
      if (is.null(sizeFactors(dds)) & is.null(normalizationFactors(dds))) {
        dds <- estimateSizeFactors(dds)
      }
      cnts <- counts(dds, normalized = normalized)
  } else if(inherits(dds, 'DESeqTransform')){
    cnts <- assay(dds)
  }

  group <- colData(dds)[[intgroup]]

  if(!all(gene %in% rownames(cnts))){
    num.to.skip <- sum(!gene %in% rownames(cnts))
    message('Skipping ', num.to.skip,
            ' genes not found in data\n')
    gene <- gene[gene %in% rownames(cnts)]
  }

  # remove '' from genes
  # NOTE: this can sometimes happen when count tables
  # are read from file, when an empty line is read into
  # the count matrix and propagates through the analysis
  gene <- setdiff(gene, '')

  # extract data for genes
  ll <- lapply(gene, function(x){
              df <- data.frame(count=cnts[x,],
                               gene=x)

              df[, intgroup] <- as.factor(group)
              df$sample <- rownames(df)
              df
          })

  ldf <- do.call(rbind, ll)

  return(ldf)
}

#' Create gene plot
#'
#' This function creates the gene plot.
#'
#' @param df data.frame with gene counts
#' @param intgroup metadata variable to plot on x-axis
#' @param factor.levels levels of intgroup to show on x-axis
#' @param title title of plot
#' @param ylab y-axis label
#' @param color metadata variable to color by
#' @param nrow number of rows to plot if faceting
#' @param ymin y-axis lower limit
#' @param ymax y-axis upper limit
#' @param log should y-axis be log10-transformed?
#' @param freey should y-axes of faceted plots have independent scales?
#' @param trendline type of trendline to draw
#' @param facet metadata variable to facet by
#' @param legend show legend?
#' @param boxes show boxes?
#' @param rotate_x_labels angle to rotate x-axis labels (default=30)
#'
#' @return ggplot handle
#'
#' @examples
#' # make example DESeq dataset
#' dds <- DESeq2::makeExampleDESeqDataSet()
#'
#' # get gene counts
#' df <- get_gene_counts(dds, gene = c('gene1', 'gene2'))
#'
#' # standard gene plot
#' p <- getcountplot(df, intgroup = "condition", factor.levels = c("A", "B"))
#'
#' # with genes faceted
#' p1 <- getcountplot(df, intgroup = "condition", factor.levels = c("A", "B"), facet = "gene")
#'
#'
#' @export
getcountplot <- function(df, intgroup='group', factor.levels, title=NULL,
                         ylab='Normalized counts', color='gene', nrow=2, ymin=NULL, ymax=NULL,
                         log=TRUE, freey=FALSE, trendline='smooth', facet=NULL, legend=TRUE, boxes=TRUE, rotate_x_labels=30){
  idx <- df[,intgroup] %in% factor.levels

  df <- df[idx,]
  df[,intgroup] <- factor(df[,intgroup], levels=factor.levels)

  # convert color and faceting variables to factor
  df[,color] <- factor(df[,color])
  if(!is.null(facet)){
    for(f in facet){
      df[,f] <- factor(df[,f])
    }
  }

  # set y-axis min & max
  ymin <- ifelse(is.null(ymin), min(df$count), ymin)
  ymax <- ifelse(is.null(ymax), max(df$count), ymax)

  p <- ggplot(df, aes(y=.data$count, x=.data[[ intgroup ]], color=.data[[ color ]], name=.data$sample)) +
    geom_point(position=position_jitterdodge(dodge.width=0.2),
                 size=2, alpha=0.5)

  if(boxes)
    p <- p + geom_boxplot(alpha=0, notch=TRUE,
                          outlier.size=0, outlier.shape=NA)

  p <- p +
    theme_bw() + ylab(ylab) + xlab('') +
    #geom_line(color='#000000') +
    theme(text = element_text(size=15),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.text.x = element_text(angle=rotate_x_labels, hjust=1))

  if(!is.null(ymax)){
    if(log){
      if(!freey){
        p <- p + scale_y_continuous(trans='log10', limits=c(ymin, ymax))
      } else {
        p <- p + scale_y_continuous(trans='log10')
      }
    } else if(!freey){
      p <- p + scale_y_continuous(limits=c(ymin, ymax))
    }
  } else if(log){
      p <- p + scale_y_continuous(trans='log10')
  }

  if(!is.null(title)) p <- p + ggtitle(title)
  if(trendline == 'smooth'){
        p <- p + geom_smooth(aes(group=.data[[ color ]]), se=FALSE, linewidth=0.5, linetype='dashed')
  } else if(trendline == 'line'){
        p <- p + stat_summary(fun=median, geom='line', linetype='dashed',
                              aes(group=.data[[ color ]]), linewidth=0.5)
  }

  if(!is.null(facet)){
    if(freey) scales <- 'free_y'
    else scales <- 'fixed'

    if(length(facet) == 1) p <- p + facet_wrap(as.formula(paste('~', facet)), nrow=nrow, scales=scales)
    else p <- p + facet_wrap(as.formula(paste('~', paste(facet, collapse=' + '))), nrow=nrow, scales=scales)
  }

  if(!legend) p <- p + theme(legend.position='none')
  p <- p + theme(axis.text.x=element_text(size=12),
                 axis.text.y=element_text(size=12),
                 axis.title.y=element_text(face='bold'))
}

#' Get initial y-axis limits
#'
#' @param df data.frame with counts. Must have column 'count'
#' @param y_delta y-axis padding for visualization, must be
#' between 0 and 1
#' @param pseudocount pseudo-count to add to the data.frame
#'
#' @return min and max limits for count column, padded for visualization
#'
#' @examples
#' # make example DESeq dataset
#' dds <- DESeq2::makeExampleDESeqDataSet()
#'
#' # get gene counts
#' df <- get_gene_counts(dds, gene = c('gene1', 'gene2'))
#'
#' # get y axis limits
#' get_y_init(df, y_delta = 0.01, pseudocount = 1)
#'
#' @export
get_y_init <- function(df, y_delta, pseudocount){
    if(!'count' %in% colnames(df))
      stop('Column "count" not found in data frame')

    if(y_delta < 0 | y_delta > 1)
      stop('y_delta must be between 0 and 1')

    df$count <- df$count + pseudocount

    df.max <- max(df$count)
    df.min <- min(df$count)

    max.init <- round((1 + y_delta)*df.max)
    min.init <- ifelse(((1 - y_delta)*df.min) > 1,
                       round((1 - y_delta)*df.min), pseudocount*0.1)

    return( c(min.init, max.init) )
}

#' Make final object for internal use by the app
#'
#' This function takes an uploaded object and sanitizes
#' it to make sure it is suitable for internal use along
#' with other additions:
#' - adds a 'dds_mapping' element that maps
#'   dds_list keys to res_list objects.
#' - if there are multiple dds_list objects, it adds
#'   a 'all_dds' element combining all samples.
#'
#' @param obj list object containing lists of DE analysis
#' results, functional enrichment objects, pattern analysis
#' objects & raw and normalized counts objects.
#'
#' @return final carnation object with additional pre-processing
#'
#' @examplesIf interactive()
#' library(DESeq2)
#'
#' # make example DESeq dataset
#' dds <- makeExampleDESeqDataSet()
#'
#' # run DE analysis
#' dds <- DESeq(dds)
#'
#' # extract comparison of interest
#' res <- results(dds, contrast = c("condition", "A", "B"))
#'
#' # perform VST normalization
#' rld <- varianceStabilizingTransformation(dds, blind = TRUE)
#'
#' # build minimal object
#' obj <- list(
#'            res_list = list(
#'                           comp = list(
#'                               res = res,
#'                               dds = "main",
#'                               label = "A vs B"
#'                           )
#'                       ),
#'            dds_list = list(main = dds),
#'            rld_list = list(main = rld)
#'        )
#'
#' # final object
#' final_obj <- make_final_object(obj)
#'
#' @export
make_final_object <- function(obj){
    # get object names and map to 'res.list', etc
    n <- names(obj)
    res.name <- n[grep('res', n)]
    dds.name <- n[grep('dds', n)]
    rld.name <- n[grep('rld', n)]
    enrich.name <- n[grep('enrich', n)]
    degpatterns.name <- n[grep('degpatterns', n)]

    # get res.list names
    comp.names <- names(obj[[res.name]])

    # if res.list contains 'res', 'dds', 'label' elements
    # build the following:
    #
    # - labels: list of labels extracted from res.list
    # - if 'symbol' column does not exist, make one from rownames
    # - dds_mapping: list mapping dds.list elements to res.list elements
    if(all(c('res','dds','label') %in% names(obj[[res.name]][[1]]))){
      labels <- lapply(obj[[res.name]], function(x){
                           x$label
          })
      names(labels) <- comp.names

      res.list <- lapply(obj[[res.name]], function(x){
        res <- x$res

        # - if no 'symbol' column, add from rownames
        # - else replace NA's with rownames
        sidx <- which('SYMBOL' %in% toupper(colnames(res)))
        if(length(sidx) == 0){
            res$symbol <- rownames(res)
        } else {
          if(any(is.na(res[,sidx]))){
              idx <- which(is.na(res[, sidx]))
              res[idx, sidx] <- rownames(res)[idx]
          }
        }
        res[order(res$padj),]
      })
      names(res.list) <- comp.names

      # make dds.list
      # latest res.list contains names of dds
      if(is.character(obj[[res.name]][[1]]$dds)){
        # keep original dds.list & rld.list
        dds.list <- obj[[dds.name]]
        rld.list <- obj[[rld.name]]

        dds_mapping <- lapply(comp.names, function(x){
                            obj[[res.name]][[x]]$dds
                        })
      } else {
        # if res.list elements contain full dds objects
        # create new dds.list & rld.list
        dds.list <- lapply(res.list, function(x) x$dds)
        names(dds.list) <- comp.names

        # make rld list
        rld.list <- lapply(comp.names, function(x){
          name <- obj[[res.name]][[x]]$dds
          obj[[rld.name]][[name]]
        })
        names(rld.list) <- comp.names

        # build dds_mapping
        dds_mapping <- as.list(comp.names)
      }
      names(dds_mapping) <- comp.names

      # replace originals with new lists
      # & add labels and dds_mapping
      obj[[res.name]] <- res.list
      obj[[dds.name]] <- dds.list
      obj[[rld.name]] <- rld.list
      obj$labels <- labels
      obj$dds_mapping <- dds_mapping

      # get symbol mapping from res objects
      all_idmap <- NULL
      all_idmap_names <- NULL

      # if symbol column exists, change rownames(dds) to symbol
      for(name in names(obj[[dds.name]])){
        # 1. try to use dds_mapping to get res.list element
        #    corresponding to dds.list object
        # 2. if not, check to see if res.list element of same
        #    name exists
        # 3. else skip
        idx <- obj$dds_mapping %in% name
        if(sum(idx) >= 1){
          res <- obj[[res.name]][[which(idx)[1]]]
        } else if(name %in% names(obj[[res.name]])){
          res <- obj[[res.name]][[name]]
        } else {
          message('no matching dds object found for ', name, ', skipping\n')
          obj[[dds.name]] <- obj[[dds.name]][!names(obj[[dds.name]]) %in% name]
          obj[[rld.name]] <- obj[[rld.name]][!names(obj[[rld.name]]) %in% name]
          next
        }

        dds <- obj[[dds.name]][[name]]
        rld <- obj[[rld.name]][[name]]

        # get symbol column from res.list object
        idx <- toupper(colnames(res)) %in% 'SYMBOL'

        # build id to symbol mapping
        if(sum(idx) > 0){
          idmap <- res[,which(idx)]

          # get gene column or use rownames
          gidx <- toupper(colnames(res)) %in% 'GENE'
          if(sum(gidx) > 0){
              idmap_names <- res[,which(gidx)]
          } else if(!is.null(rownames(res))){
              idmap_names <- rownames(res)
          } else {
              message(
                'No gene column or row names found in ',
                'res object: ', name, '. Skipping\n'
              )
              next
          }

          # change mapping to id for which symbol is NA
          idmap[is.na(idmap)] <- idmap_names[is.na(idmap)]

          # add any new genes to mapping vector
          new_idx <- which(!idmap %in% all_idmap)
          if(length(new_idx) > 0){
            all_idmap <- c(all_idmap, idmap[new_idx])
            all_idmap_names <- c(all_idmap_names, idmap_names[new_idx])
          }

          # switch rownames to symbol
          names(idmap) <- idmap_names
          rownames(dds) <- idmap[rownames(dds)]
          rownames(rld) <- idmap[rownames(rld)]
        }

        # replace NAs
        all_idmap[is.na(all_idmap)] <- all_idmap_names[is.na(all_idmap)]

        # remove duplicates from idmap
        idx <- duplicated(all_idmap)
        all_idmap <- all_idmap[which(!idx)]
        names(all_idmap) <- all_idmap_names[which(!idx)]


        obj[[dds.name]][[name]] <- dds
        obj[[rld.name]][[name]] <- rld
      }
    } else {
      # if res.list contains full DESeqDataSet objects in
      # 'dds' slot, build dds_mapping simply from res.list names
      dds_mapping <- as.list(comp.names)
      names(dds_mapping) <- comp.names

      obj$dds_mapping <- dds_mapping
    }

    # if degpatterns element exists, add symbol column
    if(length(degpatterns.name) != 0){
      obj[[ degpatterns.name ]] <- lapply(obj[[ degpatterns.name ]],
                                     function(x){
                                       # only keep 'normalized' slot
                                       if('normalized' %in% names(x)) x <- x$normalized

                                       if(!'SYMBOL' %in% toupper(colnames(x))){
                                         # if 'genes' column is found, add symbols
                                         # else do nothing
                                         gcol <- grep('GENES', toupper(colnames(x)))
                                         if(length(gcol) > 0){
                                           x$symbol <- all_idmap[x[, gcol[1]]]
                                         }
                                       }
                                       x
                                     })
    }

    # if multiple dds.list elements exist,
    # build a combined dds object to use as 'all_samples'
    #
    # Briefly, we loop over dds.list elements maintaining
    # a count matrix and a colData matrix
    # - the count matrix is subsetted to only keep genes
    #   that are common to all dds.list elements
    # - the colData matrix is subsetted to only contain
    #   columns present in all dds.list element metadata
    # - finally, the count matrix is filtered so that each
    #   sample (column) is a row in the colData matrix
    if(length(obj[[dds.name]]) > 1){
        all.counts <- NULL
        all.cdata <- NULL
        fail <- FALSE
        for(comp in names(obj[[dds.name]])){
          # skip dds.list element not present in dds_mapping
          if(!comp %in% obj$dds_mapping) next

          # initialize with first element counts & colData
          if(is.null(all.counts)){
              all.counts <- counts(obj[[dds.name]][[comp]])
              all.cdata <- colData(obj[[dds.name]][[comp]])
          } else {
              # only keep common genes (rows)
              counts.i <- counts(obj[[dds.name]][[comp]])
              common.rows <- rownames(counts.i)[rownames(counts.i) %in% rownames(all.counts)]

              # if no genes are shared, set failure flag
              if(length(common.rows) == 0){
                  fail <- TRUE
                  break
              }
              new.cols <- which(!colnames(counts.i) %in% colnames(all.counts))
              if(length(new.cols) > 0){
                  if(!all(common.rows %in% all.counts)){
                      fail <- TRUE
                      break
                  }
                  all.counts <- cbind(all.counts[common.rows,],
                                      counts.i[common.rows, new.cols])
              }

              # only keep common metadata columns
              cdata.i <- colData(obj[[dds.name]][[comp]])
              common.cols <- colnames(cdata.i)[colnames(cdata.i) %in% colnames(all.cdata)]

              # if at any point, the number of common metadata
              # columns is zero, set failure flag
              if(length(common.cols) <= 1){
                  fail <- TRUE
                  break
              }

              all.cdata <- rbind(all.cdata[, common.cols],
                                 cdata.i[colnames(counts.i)[new.cols], common.cols])
          }
        }

        # if failure flag is set (no common genes or no shared
        # metadata columns), all_dds & all_rld slots are NULL
        if(fail){
            obj$all_dds <- NULL
            obj$all_rld <- NULL
        }  else {
          # only keep sample names that are present in both
          # all.counts and all.cdata
          shared.samples <- colnames(all.counts)[colnames(all.counts) %in% rownames(all.cdata)]

          if(length(shared.samples) > 0){
            # sort all.counts cols to be same order as all.cdata rows
            all.counts <- all.counts[, shared.samples]
            all.cdata <- as.data.frame(all.cdata[rownames(all.cdata) %in% shared.samples,])

            # create dds
            all_dds <- DESeqDataSetFromMatrix(countData=as.matrix(all.counts),
                                              colData=all.cdata,
                                              design=~1)

            obj$all_dds <- all_dds
            obj$all_rld <- varianceStabilizingTransformation(all_dds, blind=TRUE)
            colData(obj$all_rld)$sample <- rownames(colData(obj$all_rld))
          } else {
            message('No shared samples in dds objects ...\n')

            obj$all_dds <- NULL
            obj$all_rld <- NULL
          }
        }
    } else if(length(obj[[dds.name]]) == 1){
        # if only one dds.list element exists, then
        # set that to all_dds & all_rld slots
        obj$all_dds <- obj[[dds.name]][[1]]
        obj$all_rld <- obj[[rld.name]][[1]]
    }

    # give final names
    all_names <- names(obj)
    new_names <- all_names
    new_names[grep(dds.name, all_names)] <- 'dds'
    new_names[grep(res.name, all_names)] <- 'res'
    new_names[grep(rld.name, all_names)] <- 'rld'
    if(length(enrich.name) > 0)
      new_names[grep(enrich.name, all_names)] <- 'enrich'
    if(length(degpatterns.name) > 0)
      new_names[grep(degpatterns.name, all_names)] <- 'degpatterns'

    names(obj) <- new_names

    return(obj)
}

#' Convert enrichResult to GeneTonic object
#'
#' This function takes an enrichResult object and
#' DE analysis results and creates a GeneTonic object.
#'
#' @param enrich enrichResult object
#' @param res data frame with DE analysis results
#'
#' @return GeneTonic object
#'
#' @examples
#' # get enrich & res objects
#' data(res_dex, package="carnation")
#' data(eres_dex, package="carnation")
#'
#' # convert to GeneTonic object
#' gt <- enrich_to_genetonic(eres_dex, res_dex)
#'
#'
#' @export
enrich_to_genetonic <- function(enrich, res){
    if(inherits(enrich, 'enrichResult'))
      l_gs <- shake_enrichResult(enrich)
    else if(inherits(enrich, 'gseaResult'))
      l_gs <- shake_gsenrichResult(enrich)

    if(!'gene' %in% colnames(res)){
      if(!is.null(rownames(res))){
        res$gene <- rownames(res)
        res <- as.data.frame(res) %>% relocate(.data$gene)
      } else {
        stop('Cannot find gene column in result data frame!')
      }
    }
    idx <- match(c('gene','symbol'), tolower(colnames(res)))
    if(length(idx) != 2){
      stop('Columns of DE results must contain "gene" & "symbol"')
    }
    anno_df <- res[,idx]
    colnames(anno_df) <- c('gene_id', 'gene_name')

    l_gs <- get_aggrscores(l_gs, res, anno_df)
    return(list(l_gs=l_gs, anno_df=anno_df))
}

#' Create a labeled MA plot
#'
#' This function creates an MA plot from a data.frame
#' containing DE analysis results.
#'
#' @param res data.frame with DE analysis results. Must contain
#'  "padj" & "log2FoldChange" columns
#' @param fdr.thres False discovery rate (FDR) threshold
#' @param fc.thres log2FoldChange threshold
#' @param fc.lim y-axis limits
#' @param lab.genes genes to label on MA plot
#' @param tolower.cols column names that will be converted to
#'  lower case
#'
#' @return ggplot handle
#'
#' @examples
#' # make mock results df
#' n_genes <- 100
#' res <- data.frame(
#'          baseMean = runif(n_genes, 10, 1000),
#'          log2FoldChange = rnorm(n_genes, 0, 2),
#'          lfcSE = runif(n_genes, 0.1, 0.5),
#'          stat = rnorm(n_genes, 0, 3),
#'          pvalue = runif(n_genes, 0, 1),
#'          padj = runif(n_genes, 0, 1),
#'          symbol = paste0("GENE", 1:n_genes),
#'          row.names = paste0("gene", 1:n_genes)
#'        )
#'
#' plotMA.label(res, lab.genes = c("gene1", "gene2"))
#'
#' @export
plotMA.label <- function(res,
                         fdr.thres=0.01,
                         fc.thres=0,
                         fc.lim=NULL,
                         lab.genes=NULL,
                         tolower.cols=c('SYMBOL','ALIAS')){
  # convert res to data frame
  res <- data.frame(res)

  if(!all(c('padj', 'log2FoldChange') %in% colnames(res))){
    stop('DE analysis results must contain "padj" & "log2FoldChange" columns')
  }

  # if y limits not specified
  if(is.null(fc.lim)){
    fc.lim <- range(res$log2FoldChange, na.rm=TRUE)
    fc.lim[1] <- floor(fc.lim[1])
    fc.lim[2] <- ceiling(fc.lim[2])
  }

  # change specific colnames to lower case
  idx <- colnames(res) %in% tolower.cols
  colnames(res)[idx] <- tolower(colnames(res)[idx])

  # change NA symbols to gene ID
  if('symbol' %in% colnames(res)){
    df <- res %>%
      mutate(geneid=rownames(res)) %>%
      mutate(symbol=as.character(.data$symbol)) %>%
      mutate(symbol = replace(.data$symbol, is.na(.data$symbol), .data$geneid[is.na(.data$symbol)]))
  } else {
    df <- res
    df$symbol <- rownames(res)
  }

  # create column with plotting character based on fc.lim
  # change plotted values for those outside plot limits
  # to the limits
  df$shape <- 'in'
  df <- df %>%
    mutate(log2FoldChange = replace(.data$log2FoldChange, is.na(.data$log2FoldChange), 0)) %>%
    mutate(log2FoldChange = replace(.data$log2FoldChange, .data$log2FoldChange > fc.lim[2], fc.lim[2])) %>%
    mutate(log2FoldChange = replace(.data$log2FoldChange, .data$log2FoldChange < fc.lim[1], fc.lim[1])) %>%
    mutate(shape = replace(.data$shape, .data$log2FoldChange == fc.lim[2], 'above')) %>%
    mutate(shape = replace(.data$shape, .data$log2FoldChange == fc.lim[1], 'below')) %>%
    mutate(shape = as.factor(.data$shape))

  # change colors for DE genes
  df$significant <- 'no'
  df <- df %>%
    mutate(significant = replace(.data$significant,
                                 .data$padj < fdr.thres & !is.na(.data$padj) & abs(.data$log2FoldChange) >= fc.thres,
                                 'yes')) %>%
    mutate(significant=as.factor(.data$significant))

  # make main plot
  p <- df %>%
    ggplot(aes(.data$baseMean, .data$log2FoldChange, color=.data$significant,
               shape=.data$shape, name=.data$symbol)) +
    geom_point(alpha=0.8) +
    ylim(fc.lim[1]-0.1, fc.lim[2]+0.1) +
    scale_x_log10()

  # add scales
  p <- p + scale_color_manual(breaks=c('no','yes'), values=c('gray40','red'), guide='legend')
  p <- p + scale_shape_manual(breaks=c('in','above','below'),
                              values=c(16, 2, 6), guide='none')

  # change theme and add horizontal line
  p <- p + theme_bw() +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
  p <- p + geom_hline(yintercept = 0, col="red", linewidth=2, alpha=0.5)


  if(!is.null(lab.genes)){
    # get data frame of genes to be labeled
    lab.list <- df %>% filter(.data$symbol %in% lab.genes)

    # add labels
    p <- p + geom_point(data=lab.list, col="black", pch=1, size=3)
    p <- p + geom_label_repel(data=lab.list, aes(label=.data$symbol),
                              fontface="italic", show.legend=FALSE,
                              max.overlaps=Inf)
  }
  return(p)

}

#' Add set column to UpSet plot matrix
#'
#' This function adds a column denoting set number
#' to a matrix generated for an upset plot with
#' fromList.with.names()
#'
#' @param df binary matrix where row = genes & columns
#'  are gene sets, with 1 indicating that a gene is present
#'  is that gene set and vice-versa
#'
#' @return data.frame with added set column
#'
#' @examples
#' # list of genes
#' lst <- list(group1 = c(a = "gene1", b = "gene2", c = "gene3", d = "gene4"),
#'             group2 = c(c = "gene3", d = "gene4"))
#'
#' # binarized matrix with group membership
#' df <- fromList.with.names(lst)
#'
#' # matrix with added set column
#' ldf <- add.set.column(df)
#'
#' @export
add.set.column <- function(df){
    # save symbol column if any
    sidx <- which(tolower(colnames(df)) %in% 'symbol')
    scol <- df[, sidx]
    df <- df %>% select(-all_of(sidx))

    # add column with set

    # first convert matrix to vector by collapsing rows, e.g.
    #
    #   0 1 0 1
    #   0 0 1 1
    #   0 0 1 1
    #
    # is converted to
    #
    #   c('0 1 0 1', '0 0 1 1', '0 0 1 1')
    df2 <- apply(df, 1, function(x) paste(x, collapse=' '))

    # count the different sequences & order decreasing
    set.counts <- table(df2)
    set.counts <- set.counts[order(set.counts, decreasing=TRUE)]

    # map sequences to set names, e.g. set1, set2
    #
    # NOTE:
    # - set numbers are padded with 0 to have the same number of digits
    #   e.g. if there are 10 sets, sets are numbered 'set01' ... 'set10'
    max.digits <- nchar(length(set.counts))
    set.names <- unlist(lapply(seq_len(length(set.counts)),
                               function(x){
                                   paste0('set', paste(rep(0, max.digits-nchar(x)), collapse=''), x)
                               }))
    names(set.names) <- names(set.counts)

    # now add 'set' column with genes mapped to set name
    # then order by set
    df$set <- set.names[df2]
    df$symbol <- scol

    df <- df[order(df$set),]

    # move set & symbol column to beginning
    df <- df %>% relocate('set') %>%
      relocate('symbol')

    return(df)
}

#' Create an interactive labeled MA plot
#'
#' This function creates an MA plot from a data.frame
#' containing DE analysis results using plot_ly
#'
#' @param res data.frame with DE analysis results. Must contain
#'  "padj" & "log2FoldChange" columns
#' @param fdr.thres False discovery rate (FDR) threshold
#' @param fc.thres log2FoldChange threshold
#' @param fc.lim y-axis limits
#' @param lab.genes genes to label on MA plot
#' @param tolower.cols column names that will be converted to
#'  lower case
#'
#' @return plotly handle
#'
#' @examples
#' # make mock results df
#' n_genes <- 100
#' res <- data.frame(
#'          baseMean = runif(n_genes, 10, 1000),
#'          log2FoldChange = rnorm(n_genes, 0, 2),
#'          lfcSE = runif(n_genes, 0.1, 0.5),
#'          stat = rnorm(n_genes, 0, 3),
#'          pvalue = runif(n_genes, 0, 1),
#'          padj = runif(n_genes, 0, 1),
#'          symbol = paste0("GENE", 1:n_genes),
#'          row.names = paste0("gene", 1:n_genes)
#'        )
#'
#' plotMA.label_ly(res, lab.genes = c("gene1", "gene2"))
#'
#' @export
plotMA.label_ly <- function(res,
                         fdr.thres=0.01,
                         fc.thres=0,
                         fc.lim=NULL,
                         lab.genes=NULL,
                         tolower.cols=c('SYMBOL','ALIAS')){
  # convert res to data frame
  res <- data.frame(res)

  if(!all(c('padj', 'log2FoldChange') %in% colnames(res))){
    stop('DE analysis results must contain "padj" & "log2FoldChange" columns')
  }

  # if y limits not specified
  if(is.null(fc.lim)){
    fc.lim <- range(res$log2FoldChange, na.rm=TRUE)
    fc.lim[1] <- floor(fc.lim[1])
    fc.lim[2] <- ceiling(fc.lim[2])
  }

  # change specific colnames to lower case
  idx <- colnames(res) %in% tolower.cols
  colnames(res)[idx] <- tolower(colnames(res)[idx])

  # change NA symbols to gene ID
  if('symbol' %in% colnames(res)){
    df <- res %>%
      mutate(geneid=rownames(res)) %>%
      mutate(symbol=as.character(.data$symbol)) %>%
      mutate(symbol = replace(.data$symbol, is.na(.data$symbol), .data$geneid[is.na(.data$symbol)]))
  } else {
    df <- res
    df$symbol <- rownames(res)
  }

  # create column with plotting character based on fc.lim
  # change plotted values for those outside plot limits
  # to the limits
  df$shape <- 'in'
  df <- df %>%
    mutate(log2FoldChange = replace(.data$log2FoldChange, is.na(.data$log2FoldChange), 0)) %>%
    mutate(log2FoldChange = replace(.data$log2FoldChange, .data$log2FoldChange > fc.lim[2], fc.lim[2])) %>%
    mutate(log2FoldChange = replace(.data$log2FoldChange, .data$log2FoldChange < fc.lim[1], fc.lim[1])) %>%
    mutate(shape = replace(.data$shape, .data$log2FoldChange == fc.lim[2], 'above')) %>%
    mutate(shape = replace(.data$shape, .data$log2FoldChange == fc.lim[1], 'below')) %>%
    mutate(shape = as.factor(.data$shape))

  # change colors for DE genes
  df$significant <- 'no'
  df <- df %>%
    mutate(significant = replace(.data$significant,
                                 .data$padj < fdr.thres & !is.na(.data$padj) & abs(.data$log2FoldChange) >= fc.thres,
                                 'yes')) %>%
    mutate(significant=as.factor(.data$significant))

  # get n.s. points first
  ns.df <- df %>% filter(.data$significant == 'no')
  ns.below <- ns.df %>%
      filter(.data$log2FoldChange == fc.lim[1])
  ns.above <- ns.df %>%
      filter(.data$log2FoldChange == fc.lim[2])
  ns.rest <- ns.df %>%
      filter(.data$log2FoldChange > fc.lim[1] & .data$log2FoldChange < fc.lim[2])

  # next get de points
  de.df <- df %>% filter(.data$significant == 'yes')
  de.below <- de.df %>%
      filter(.data$log2FoldChange == fc.lim[1])
  de.above <- de.df %>%
      filter(.data$log2FoldChange == fc.lim[2])
  de.rest <- de.df %>%
      filter(.data$log2FoldChange > fc.lim[1] & .data$log2FoldChange < fc.lim[2])

  # make main plot with plot_ly
  #
  # The main idea here is to make separate layers for
  # each class of points & remove hoverinfo for n.s.
  # points.
  p <- plot_ly(x=ns.rest$baseMean, y=ns.rest$log2FoldChange, type='scatter',
              text=ns.rest$symbol,
              mode='markers',
              hoverinfo='text',
              name='no',
              marker=list(color='gray', alpha=0.3)) %>%
        layout(xaxis=list(type='log',
                          title='baseMean',
                          showgrid=FALSE),
               yaxis=list(title='log2FoldChange',
                          showgrid=FALSE,
                          range=c(fc.lim[1]-0.1, fc.lim[2]+0.1)),
               legend=list(title=list(text='<b> significant </b>')))

  # add not DE points below y limits
  if(nrow(ns.below) > 0){
   p <- p %>% add_markers(x=ns.below$baseMean, y=ns.below$log2FoldChange,
                  text=ns.below$symbol,
                  hoverinfo='text',
                  marker=list(
                      color='gray', size=5,
                      alpha=0.3,
                      symbol='triangle-down-open'),
                  showlegend=FALSE)
  }

  # add not DE points above y limits
  if(nrow(ns.above) > 0){
   p <- p %>% add_markers(x=ns.above$baseMean, y=ns.above$log2FoldChange,
                  text=ns.above$symbol,
                  hoverinfo='text',
                  marker=list(
                      color='gray', size=5,
                      alpha=0.3,
                      symbol='triangle-up-open'),
                  showlegend=FALSE)
  }

  # add DE points below y limits
  if(nrow(de.below) > 0){
   p <- p %>% add_markers(x=de.below$baseMean, y=de.below$log2FoldChange,
                  text=de.below$symbol,
                  hoverinfo='text',
                  marker=list(
                      color='red', size=5,
                      alpha=0.3,
                      symbol='triangle-down-open'),
                  showlegend=FALSE)
  }

  # add DE points above y limits
  if(nrow(de.above) > 0){
   p <- p %>% add_markers(x=de.above$baseMean, y=de.above$log2FoldChange,
                  text=de.above$symbol,
                  hoverinfo='text',
                  marker=list(
                      color='red', size=5,
                      alpha=0.3,
                      symbol='triangle-up-open'),
                  showlegend=FALSE)
  }

  # add DE points within y limits
  if(nrow(de.rest) > 0){
    p <- p %>%
        add_markers(x=de.rest$baseMean, y=de.rest$log2FoldChange,
                    name='yes',
                    text=de.rest$symbol,
                    hoverinfo='text',
                    marker=list(color='red', alpha=0.3))
  }


  # add labels if any
  if(!is.null(lab.genes)){
    # get data frame of genes to be labeled
    lab.list <- df %>% filter(.data$symbol %in% lab.genes)

    if(nrow(lab.list) > 0){
      p <- p %>%
          add_markers(x=lab.list$baseMean, y=lab.list$log2FoldChange,
                      text=lab.list$symbol,
                      hoverinfo='none',
                      name='labeled',
                      marker=list(color='black',
                                  symbol='circle-open',
                                  size=10,
                                  line=list(width=2)))
    }
  }

  # add horizontal y = 0 line
  p <- p %>%
      layout(shapes=list(type = "line",
                         x0 = 0, x1 = 1,
                         xref = "paper",
                         y0 = 0, y1 = 0,
                         line = list(color = 'red', width=3, alpha=0.3)))

  return(p)

}

#' Get top DE genes by log2FoldChange or adjusted p-value
#'
#' @param res data.frame with DE analysis results
#' @param fdr.thres FDR threshold
#' @param fc.thres log2FoldChange threshold
#' @param n number of genes to return
#' @param by metric to determine top genes ('log2FoldChange' or 'padj')
#'
#' @return vector of gene symbols
#'
#' @examples
#' # get DE results
#' data(res_dex, package='carnation')
#'
#' g <- top.genes(res_dex)
#
#' @export
top.genes <- function(res, fdr.thres=0.01, fc.thres=0, n=10, by='log2FoldChange'){

  # sort by padj
  res <- res[order(res$padj),]

  # DE genes
  idx <- res$padj < fdr.thres & !is.na(res$padj) & abs(res$log2FoldChange) >= fc.thres
  if(sum(idx, na.rm=TRUE) == 0) return( NULL )
  res.de <- res[idx,]

  # get symbol column, if any
  if('symbol' %in% tolower(colnames(res))){
    sidx <- which(tolower(colnames(res)) %in% c('symbol'))

    # remove NAs if any
    na.idx <- is.na(res.de[, sidx])
    res.de[which(na.idx), sidx] <- rownames(res.de)[which(na.idx)]
  } else {
    sidx <- NULL
  }

  # number of upregulated genes
  if(nrow(res.de) == 0) return( NULL )
  else if(nrow(res.de) <= n) n <- nrow(res.de)

  if(by == 'log2FoldChange'){
    # order DE genes by lfc
    res.de <- res.de[order(res.de$log2FoldChange, decreasing=TRUE),]

    n.up <- ceiling(n/2)

    # get top upregulated & downregulated genes
    top.idx <- c(seq_len(n.up), (nrow(res.de) - (n - n.up) + 1):nrow(res.de))
  } else if(by == 'padj'){
    top.idx <- seq_len(n)
  }

  if(is.null(sidx)) return( unique(rownames(res.de)[top.idx]) )
  else return( unique(res.de[top.idx, sidx]) )
}

#' Add metadata to counts data frame
#'
#' @param df data.frame with gene counts
#' @param coldata data.frame with metadata
#' @param exclude.intgroups metadata columns to ignore
#'
#' @return counts data frame with added metadata
#'
#' @examples
#' library(DESeq2)
#'
#' # make example DESeq data set
#' dds <- makeExampleDESeqDataSet()
#'
#' # extract counts and metadata
#' df <- assay(dds)
#' coldata <- colData(dds)
#'
#' # get gene counts df
#' counts_df <- get_gene_counts(dds, paste0('gene', seq_len(10)))
#'
#' # add metadata
#' counts_df <- add_metadata(counts_df, coldata, exclude.intgroups=NULL)
#'
#' @export
add_metadata <- function(df, coldata, exclude.intgroups){
  coldata$sample <- rownames(coldata)

  coldata <- coldata %>% as.data.frame %>% select(-any_of(exclude.intgroups))

  if(any(colnames(df) %in% colnames(coldata))){
    common <- intersect(colnames(df), colnames(coldata))

    df <- merge(df, coldata, by=common)
  } else {
    df <- merge(df, coldata, by='sample')
  }

  return(df)
}

#' format gene names to look pretty in table output
#'
#' This function works by grouping long lists of genes
#' into groups of a specified size. Each group is collapsed
#' using commas, while groups are separated by spaces
#' so that datatable formatting is tricked into separating
#' space-separated groups and not comma-separated groups
#'
#' @param g vector of gene names
#' @param sep gene name separator
#' @param genes.per.line number of genes to show in a line
#'
#' @return vector of gene names prettified for data.table output
#'
#' @examples
#' # string with genes separated by '/'
#' g <- "gene1/gene2/gene3/gene4/gene5/gene6/gene7"
#'
#' gg <- format_genes(g, genes.per.line=3)
#'
#' @export
format_genes <- function(g, sep='\\/', genes.per.line=6){
  # get number of genes associated with each functional term
  g2 <- lapply(g, function(x) unlist(strsplit(x, split=sep)))
  ngenes <- unlist(lapply(g2, length))

  # get idx of terms with ngenes > genes.per.line
  idx <- ngenes >= genes.per.line

  # if no hits return prev vector with sep replaced
  if(length(idx) == 0) return(unlist(gsub(sep,',',g)))
  else {
    final <- g
    # for terms with ngenes <= genes.per.line,
    # just collapse with new sep
    final[which(!idx)] <- unlist(gsub(sep,',',g[which(!idx)]))

    final[which(idx)] <- unlist(lapply(which(idx),
                    function(x){
                      g.tmp <- g2[[x]]
                      nsets <- floor(length(g.tmp)/genes.per.line)
                      rem <- length(g.tmp) %% genes.per.line

                      if(rem == 0){
                        start.idx <- (0:(nsets-1))*genes.per.line + 1
                        end.idx <- (seq_len(nsets))*genes.per.line
                      } else {
                        start.idx <- (0:nsets)*genes.per.line + 1
                        end.idx <- c((seq_len(nsets))*genes.per.line, length(g.tmp))
                      }
                      gg <- unlist(lapply(seq_len(length(start.idx)),
                                function(x){
                                    tmp <- paste(g.tmp[start.idx[x]:end.idx[x]],
                                                 collapse=',')
                                    paste0(tmp, ',')
                                }))
                      paste(gg, collapse=' ')
                    }))
  }
  final
}

#' Prepare list for UpSet plots, but include rownames
#'
#' @param lst List of sets to compare (same input as to UpSetR::fromList)
#'
#' @return data.frame of 1 and 0 showing which genes are in which sets
#'
#' @examples
#' # list of genes
#' lst <- list(group1 = c(a = "gene1", b = "gene2", c = "gene3", d = "gene4"),
#'             group2 = c(c = "gene3", d = "gene4"))
#'
#' # binarized matrix with group membership
#' df <- fromList.with.names(lst)
#'
#' @export
fromList.with.names <- function(lst){
  element_names <- unique(do.call('rbind',
                     lapply(lst, function(x){
                       data.frame(id=names(x),
                                  symbol=unname(x),
                                  row.names=NULL)
                     })
                   ))

  data <- unlist(lapply(lst, function(x) {
    x <- as.vector(match(element_names$id, names(x)))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(lst), byrow = FALSE))

  # This is the only way in which this function differs from UpSetR::fromList
  # NOTE: here we use the unique column as rownames
  if(sum(duplicated(element_names$id)) == 0)
    rownames(data) <- element_names$id
  else
    rownames(data) <- element_names$symbol

  data <- data[which(rowSums(data) != 0), ]
  colnames(data) <- names(lst)

  # add symbol column
  # NOTE: if column 'id' was not unique, then 'symbol' is used instead
  if(sum(duplicated(element_names$id)) == 0)
    data$symbol <- element_names$symbol
  else
    data$symbol <- element_names$id
  data <- data %>% relocate(.data$symbol)

  return(data)
}

#' Plot a degPatterns object
#'
#' This function plots a degPatterns object.
#'
#' @param obj degPatterns object
#' @param time metadata variable to plot on x-axis
#' @param color variable to color plot
#' @param cluster_column column to use for grouping genes
#' @param cluster_to_show which clusters to show in plot
#' @param x_order order of x-axis values
#' @param points boolean, show samples on plot? Default: TRUE
#' @param boxes boolean, show boxes on plot? Default: TRUE
#' @param smooth what type of trendline to use? can be 'smooth' (default) or 'line'.
#' @param lines show lines joining samples? Default: TRUE
#' @param facet boolean, should plot be faceted? Default: TRUE
#' @param prefix_title string, prefix for facet titles
#' @param genes_to_label genes to label on plot
#'
#' @return ggplot handle
#'
#' @examples
#' # get degpatterns object
#' data(degpatterns_dex, package = 'carnation')
#'
#' # get pattern plot
#' all_clusters <- unique(degpatterns_dex$normalized$cluster)
#'
#' dp <- get_degplot(degpatterns_dex, time='dex',
#'                   cluster_to_show=all_clusters,
#'                   x_order=c('untrt','trt'))
#'
#' @export
get_degplot <- function(obj, time, color=NULL,
                        cluster_column='cluster',
                        cluster_to_show,
                        x_order,
                        points=TRUE,
                        boxes=TRUE,
                        smooth='smooth',
                        lines=TRUE,
                        facet=TRUE,
                        prefix_title='Cluster ',
                        genes_to_label=NULL){
  # normalized element from degpatterns output or dataframe
  if(is.data.frame(obj)) table <- obj
  else table <- obj$normalized

  # subset n based on clusters of size > minc
  nodup.idx <- !duplicated(table$genes)
  cluster.sizes <- table(table[nodup.idx, cluster_column])
  cluster.sizes <- cluster.sizes[names(cluster.sizes) %in% cluster_to_show]

  # this part taken from degPlotCluster & tweaked
  if (cluster_column  %in% colnames(table)){
      table[['cluster']] <- table[[cluster_column]]
      counts <- cluster.sizes

      table <- inner_join(table,
                 data.frame(cluster = as.integer(names(counts)),
                            title = paste0(prefix_title,
                                           names(counts),
                                           " (n = " ,
                                           counts,
                                           ")"),
                            stringsAsFactors = FALSE),
                 by = "cluster")
  }

  if (is.null(color)){
      color <- "dummy"
      table[[color]] <- "one_group"
  }

  if("symbol" %in% colnames(table)){
    table[["line_group"]] <- paste(table[["symbol"]],
                                table[[color]])
  } else {
    table[["line_group"]] <- paste(table[["genes"]],
                                  table[[color]])
  }

  # get time variable in specific order
  idx <- table[[time]] %in% x_order

  table <- table[idx,]
  table[[time]] <- factor(table[[time]], levels=x_order)

  if(!is.null(genes_to_label)){
    if('symbol' %in% colnames(table)){
        gidx <- table[['symbol']] %in% genes_to_label

        # subset to keep labeled genes, then unique to count
        label_table <- table[gidx,]
        label_table <- label_table[!duplicated(label_table[, 'symbol']),]
    } else {
        gidx <- table[['genes']] %in% genes_to_label

        # subset to keep labeled genes, then unique to count
        label_table <- table[gidx,]
        label_table <- label_table[!duplicated(label_table[, 'genes']),]
    }

    # count number of labeled genes in clusters,
    # then attach to title for that cluster
    label_counts <- table(label_table[, 'title'])
    for(cl in names(label_counts)){
      idx <- table[ ,'title'] == cl
      table[idx, 'title'] <- paste0(table[idx, 'title'], ', ',
                                    label_counts[cl], ' labeled')
    }

    # - first get the cluster -> title mapping
    # - then, order title by cluster order
    # NOTE: doing this here since labels change titles
    title_df <- unique(table[, c('cluster', 'title')])
    cl_titles <- title_df[[ 'title' ]]
    names(cl_titles) <- title_df[[ 'cluster' ]]
    table[[ 'title' ]] <- factor(table[['title']], levels=cl_titles[ cluster_to_show ])

    # do this again to update values of title
    label_table <- table[gidx,]

    # if labeling genes background lines become grayed out
    p <- ggplot(table, aes(x = .data[[ time ]], y = .data$value,
                           fill = .data[[ color ]]))


    # don't show lines/points if labeling genes
    lines <- FALSE
    points <- FALSE
    base_color <- 'lightgray'

    if (boxes)
        p <- p + geom_boxplot(alpha = 0,
                              #aes_string(fill=color),
                              color = base_color,
                              outlier.size = 0,
                              outlier.shape = NA)

    splan <- length(unique(table[[time]])) - 1L

    if(smooth == 'smooth'){
      p <- p + geom_smooth(aes(x=.data[[ time ]], y=.data$value,
                               group=.data[[ color ]]),
                           color=base_color,
                           se=FALSE)
                           #method = "lm",
                           #formula = y~poly(x, splan))
    } else if(smooth == 'line'){
      p <- p + stat_summary(
                  fun=median,
                  geom='line',
                  aes(x=.data[[ time ]], y=.data$value,
                      group=.data[[ color ]]),
                  color=base_color,
                  linewidth=1)
    }

    if (facet)
        p <- p + facet_wrap(~title)

    # turn off fill guide
    p <- p + guides(fill='none')

    # check to make sure that some labeled genes made it
    # through after filtering
    if(length(unique(label_table$line_group)) > 10){
      if(smooth == 'smooth'){
        p <- p + geom_smooth(data=label_table,
                           aes(group=.data$line_group,
                               linetype=color),
                           color='#F8766D',
                           alpha=0.5,
                           se=FALSE,
                           linewidth=0.75)
      } else if(smooth == 'line'){
        p <- p + stat_summary(data=label_table,
                           fun=median,
                           geom='line',
                           aes(group=.data$line_group,
                               linetype=color),
                           color='#F8766D',
                           alpha=0.5,
                           linewidth=0.75)
      }
    } else if(nrow(label_table) > 0){
      if(smooth == 'smooth'){
        p <- p + geom_smooth(data=label_table,
                           aes(group=.data$line_group,
                               color=.data$symbol,
                               linetype=.data[[ color ]]),
                           se=FALSE,
                           linewidth=0.75)
      } else if(smooth == 'line'){
        p <- p + stat_summary(data=label_table,
                           fun=median,
                           geom='line',
                           aes(group=.data$line_group,
                               color=.data$symbol,
                               linetype=.data[[ color ]]),
                           linewidth=0.75)
      }
    }
  } else {
    # - first get the cluster -> title mapping
    # - then, order title by cluster order
    title_df <- unique(table[, c('cluster', 'title')])
    cl_titles <- title_df[[ 'title' ]]
    names(cl_titles) <- title_df[[ 'cluster' ]]
    table[[ 'title' ]] <- factor(table[['title']], levels=cl_titles[ cluster_to_show ])

    p <- ggplot(table, aes(x = .data[[ time ]], y = .data$value,
                           fill = .data[[ color ]],
                           color = .data[[ color ]]))

    if (boxes)
        p <- p + geom_boxplot(alpha = 0,
                              outlier.size = 0,
                              outlier.shape = NA)
    if (points)
        p <- p +
        geom_point(alpha = 0.4, size = 1,
                   position = position_jitterdodge(dodge.width = 0.9))

    splan <- length(unique(table[[time]])) - 1L

    if(smooth == 'smooth'){
          p <- p + geom_smooth(aes(x=.data[[ time ]], y=.data$value,
                                   color=.data[[ color ]],
                                   group=.data[[ color ]]),
                               se=FALSE)
                               #method = "lm",
                               #formula = y~poly(x, splan))
    } else if(smooth == 'line'){
      p <- p + stat_summary(
                  fun=median,
                  geom='line',
                  aes(x=.data[[ time ]], y=.data$value,
                      color=.data[[ color ]],
                      group=.data[[ color ]]),
                  linewidth=1)
    }

    if (lines){
        p <- p + geom_line(aes(group = .data$line_group),
                               alpha = 0.1)
    }

    if (facet)
        p <- p + facet_wrap(~title)

  }

  p <- p +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90,
                                       hjust = 1,
                                       vjust = 0.5),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            text=element_text(size=15)) +
      ylab("Z-score of gene abundance") +
      xlab("")

  return(p)
}

#' Combine everything in the results list into a single table
#'
#' @param res.list Named list of lists, where each sublist contains the following
#'                 names: c('res', 'dds', 'label'). "res" is a DESeqResults object,
#'                 "dds" is either the indexing label for the dds.list object or
#'                  the DESeq object, and "label" is a nicer-looking
#'                 label to use. NOTE: backwards compatibility with older versions
#'                  of lcdb-wf depends on no dds.list object being passed.
#' @param dds.list List of DESeqDataSet objects whose names are
#'                 expected to match 'dds' slots in the 'res.list'
#'                 object
#' @param dds_mapping List mapping names of dds.list to res.list
#'                    elements
#' @param alpha false-discovery rate threshold
#' @param lfc.thresh log2FoldChange threshold
#' @param labels list of descriptions for res.list elements
#'
#' @return Dataframe
#'
#' @examples
#' n_genes <- 100
#'
#' #  make mock dds list
#' dds_list <- list(main=DESeq2::makeExampleDESeqDataSet(n=n_genes))
#'
#' # make mock results df
#' res1 <- data.frame(
#'           baseMean = runif(n_genes, 10, 1000),
#'           log2FoldChange = rnorm(n_genes, 0, 2),
#'           lfcSE = runif(n_genes, 0.1, 0.5),
#'           stat = rnorm(n_genes, 0, 3),
#'           pvalue = runif(n_genes, 0, 1),
#'           padj = runif(n_genes, 0, 1),
#'           symbol = paste0("GENE", 1:n_genes),
#'           row.names = paste0("gene", 1:n_genes)
#'         )
#'
#' res2 <- data.frame(
#'           baseMean = runif(n_genes, 10, 1000),
#'           log2FoldChange = rnorm(n_genes, 0, 2),
#'           lfcSE = runif(n_genes, 0.1, 0.5),
#'           stat = rnorm(n_genes, 0, 3),
#'           pvalue = runif(n_genes, 0, 1),
#'           padj = runif(n_genes, 0, 1),
#'           symbol = paste0("GENE", 1:n_genes),
#'           row.names = paste0("gene", 1:n_genes)
#'         )
#'
#' # make list of results
#' res_list <- list(
#'               comp1=res1,
#'               comp2=res2
#'             )
#'
#' # make dds mapping
#' dds_mapping <- list(comp1='main', comp2='main')
#'
#' # get summary
#' df <- summarize.res.list(res_list, dds_list, dds_mapping, alpha=0.1, lfc.thresh=0)
#'
#' @export
summarize.res.list <- function(res.list, dds.list, dds_mapping, alpha, lfc.thresh, labels=NULL){
    slist <- list()
    for (name in names(res.list)){
        x <- my.summary(res.list[[name]], dds.list[[ dds_mapping[[name]] ]], alpha, lfc.thresh)
        if(!is.null(labels)){
            slist[[name]] <- cbind('comparison'=name, 'description'=labels[[name]], x)
        } else {
            slist[[name]] <- cbind('comparison'=name, x)
        }
    }
    slist <- do.call(rbind, slist)
    rownames(slist) <- names(res.list)
    return(slist)
}

#' Summarize DESeq2 results into a dataframe
#'
#' summary(res) prints out info; this function captures it into a dataframe
#'
#' @param res DESeq2 results object
#' @param dds DEseq2 object
#' @param alpha Alpha level at which to call significantly changing genes
#' @param lfc.thresh log2FoldChange threshold
#'
#' @return Dataframe of summarized results
#'
#' @examples
#' n_genes <- 100
#'
#' #  make mock dds list
#' dds <- DESeq2::makeExampleDESeqDataSet(n=n_genes)
#'
#' # make mock results df
#' res <- data.frame(
#'          baseMean = runif(n_genes, 10, 1000),
#'          log2FoldChange = rnorm(n_genes, 0, 2),
#'          lfcSE = runif(n_genes, 0.1, 0.5),
#'          stat = rnorm(n_genes, 0, 3),
#'          pvalue = runif(n_genes, 0, 1),
#'          padj = runif(n_genes, 0, 1),
#'          symbol = paste0("GENE", 1:n_genes),
#'          row.names = paste0("gene", 1:n_genes)
#'        )
#'
#' # get summary
#' df <- my.summary(res, dds, alpha=0.1)
#'
#' @export
my.summary <- function(res, dds, alpha, lfc.thresh=0){
   #if (missing(alpha)){
   #    alpha <- if (is.null(metadata(res)$alpha)){ 0.1 } else { metadata(res)$alpha }
   #}
   notallzero <- sum(res$baseMean > 0)
   up <- sum(res$padj < alpha & res$log2FoldChange > lfc.thresh, na.rm=TRUE)
   down <- sum(res$padj < alpha & res$log2FoldChange < -lfc.thresh, na.rm=TRUE)
   filt <- sum(!is.na(res$pvalue) & is.na(res$padj))
   outlier <- sum(res$baseMean > 0 & is.na(res$pvalue))
   #ft <- if(is.null(metadata(res)$filterThreshold)){ 0 } else { round(metadata(res)$filterThreshold) }

   if(inherits(res, 'DESeqResults')){
     # figure out contrast from res@elementMetadata$description
     lfc.idx <- grep('log2 fold change', res@elementMetadata$description)
     contrast <- res@elementMetadata$description[lfc.idx]

     # parse contrast string for prettier display
     # remove any character ('.+') before colon ('\\:') and trailing space (\\s)
     contrast <- sub('.+\\:\\s*', '', contrast)
   } else {
     contrast <- NA
   }

   if(!is.null(dds)) curr_design <- deparse(design(dds), width.cutoff=500L)
   else curr_design <- NA

   # adjust width.cutoff as newline insertion causes this to return a df with
   # multiple duplicate rows!
   df <- data.frame(
                    #alpha=alpha,
                    #lfcThreshold=lfc.thresh,
                    up=up,
                    down=down,
                    total.genes=nrow(res),
                    total.nonzero=notallzero,
                    outliers=outlier,
                    low.counts=filt,
                    design=curr_design,
                    contrast=contrast
                    )
   return(df)
}

#' Plot an interactive PCA plot
#'
#' @param rld DESeqTransform object output by varianceStabilizingTransformation() or rlog()
#' @param intgroup character vector of names in colData(x) to use for grouping
#'
#' @return Handle to ggplot with added label field in aes_string() for plotting with ggplotly()
#'
#' @examples
#' # make example dds object
#' dds <- DESeq2::makeExampleDESeqDataSet()
#'
#' # normalize
#' rld <- DESeq2::varianceStabilizingTransformation(dds, blind=TRUE)
#'
#' # make pca plot
#' p <- plotPCA.ly(rld, intgroup='condition')
#'
#' @export
plotPCA.ly <- function(rld, intgroup){
  mat <- plotPCA(rld, intgroup, returnData=TRUE)
  pv <- attr(mat, 'percentVar')
  p <- ggplot(data=mat, aes(x=.data$PC1, y=.data$PC2, color=.data$group, label=.data$name)) +
    geom_point(size=3) + xlab(paste0('PC1: ', round(pv[1]*100), '% variance')) +
    ylab(paste0('PC2: ', round(pv[2]*100), '% variance')) + coord_fixed()
  return(p)
}

############### GENETONIC FUNCTIONS #######################

#' Radar plot
#'
#' This is a copy of gs_radar from GeneTonic where the labels of
#' gene sets are converted to parameters
#'
#' @param res_enrich GeneTonic object for comparison 1
#' @param res_enrich2 GeneTonic object for comparison 2 (default = NULL)
#' @param label1 label for comparison 1
#' @param label2 label for comparison 2
#' @param n_gs number of gene sets (default = 20)
#' @param p_value_column column to use as p-value (default = 'gs_pvalue')
#'
#' @return ggplot handle
#'
#' @examples
#' library(GeneTonic)
#'
#' # get DESeqResults object
#' data(res_dex, package='carnation')
#'
#' # get enrichResult object
#' data(eres_dex, package='carnation')
#'
#' # convert to GeneTonic object
#' gt <- shake_enrichResult(eres_dex)
#'
#' # get annotation df
#' idx <- match(c('gene','symbol'), tolower(colnames(res_dex)))
#' anno_df <- res_dex[,idx]
#' colnames(anno_df) <- c('gene_id', 'gene_name')
#'
#' # add aggregate score columns
#' gt <- get_aggrscores(gt, res_dex, anno_df)
#'
#' # make radar plot
#' p <- gs_radar(gt)
#'
#' @export
gs_radar <- function(res_enrich,
                     res_enrich2 = NULL,
                     label1 = 'scenario 1',
                     label2 = 'scenario 2',
                     n_gs = 20,
                     p_value_column = "gs_pvalue") {

  # res_enrich has to contain the Z-score to be displayed
  if (!("z_score" %in% colnames(res_enrich))) {
    warning("You need to add the z_score or the aggregated score")
  }

  if (!is.null(res_enrich2)) {
    if (!("z_score" %in% colnames(res_enrich2))) {
      warning("You need to add the z_score or the aggregated score")
    }
  }

  # only one set
  if (is.null(res_enrich2)) {
    res_enrich$logp10 <- -log10(res_enrich[[p_value_column]])
    res_enrich <- res_enrich[seq_len(n_gs), ]
    log_smallest_p <- max(res_enrich$logp10)
    set_colors <- brewer.pal(n = 8, "Set1")

    p <- plot_ly(
      type = "scatterpolar",
      mode = "markers",
      fill = "toself"
    ) %>%
      add_trace(
        r = c(res_enrich$logp10, res_enrich$logp10[1]), # recycling the first element
        theta = c(res_enrich[["gs_description"]], res_enrich[["gs_description"]][1]),
        name = label1
      ) %>%
      layout(
        polar = list(radialaxis = list(visible = TRUE,
                                       range = c(0, log_smallest_p))
        )
        # ,
        # title = "Geneset Radar Chart", font = list(size = 10)
      )
  } else {
    # if res_enrich2 is also provided
    gs_set1 <- res_enrich$gs_id
    gs_set2 <- res_enrich2$gs_id
    gs_common <- intersect(gs_set1, gs_set2)
    # restrict to the top common n_gs
    gs_common <- gs_common[seq_len(min(n_gs, length(gs_common)))]

    if (length(gs_common) == 0) {
      stop("No gene sets have been found in common to the two enrichment results")
    }

    common_re1 <- res_enrich[gs_common, ]
    common_re2 <- res_enrich2[gs_common, ]

    common_re1$logp10 <- -log10(common_re1[[p_value_column]])
    common_re2$logp10 <- -log10(common_re2[[p_value_column]])
    # if needed, I could access Z and aggregated scores

    log_smallest_p <- max(common_re1$logp10, common_re2$logp10)
    set_colors <- brewer.pal(n = 8, "Set1")

    p <- plot_ly(
      type = "scatterpolar",
      mode = "markers",
      fill = "toself"
    ) %>%
      add_trace(
        r = c(common_re1$logp10, common_re1$logp10[1]), # recycling the first element
        theta = c(common_re1[["gs_description"]], common_re1[["gs_description"]][1]),
        name = label1
      ) %>%
      add_trace(
        r = c(common_re2$logp10, common_re2$logp10[1]),
        theta = c(common_re2[["gs_description"]], common_re2[["gs_description"]][1]),
        name = label2
      ) %>%
      layout(
        polar = list(radialaxis = list(visible = TRUE,
                                       range = c(0, log_smallest_p))
        )
      )
  }

  return(p)
}

#' Adjustable PCA plot
#'
#' Create a PCA plot with specified PCs on x- and y-axis
#'
#' @param object normalized DESeqDataSet object
#' @param intgroup metadata variable to use for grouping samples
#' @param pcx principal component to plot on x-axis
#' @param pcy principal component to plot on y-axis
#' @param pcz principal component to plot on z-axis. If not NULL,
#'        function returns a 3-D PCA plot.
#' @param ntop number of most-variable genes to use
#' @param samples vector of sample names to show on plot
#' @param loadings boolean, show gene loadings? Default is FALSE.
#' @param loadings_ngenes integer, # genes to show loadings for (default=10)
#'
#' @return ggplot handle
#'
#' @examples
#' # make example dds object
#' dds <- DESeq2::makeExampleDESeqDataSet()
#'
#' # normalize
#' rld <- DESeq2::varianceStabilizingTransformation(dds, blind=TRUE)
#'
#' # make pca plot
#' p <- plotPCA.san(rld, intgroup='condition', pcx='PC1', pcy='PC2')
#'
#'
#' @export
plotPCA.san <- function (object, intgroup = "group",
                         pcx, pcy, pcz=NULL,
                         ntop = 500,
                         samples=NULL,
                         loadings=FALSE,
                         loadings_ngenes=10){
  pcx <- as.numeric(pcx)
  pcy <- as.numeric(pcy)
  if(!is.null(pcz)) pcz <- as.numeric(pcz)

  mat <- assay(object)
  cdata <- colData(object)

  if(!is.null(samples)){
    mat <- mat[, colnames(mat) %in% samples]
    cdata <- cdata[rownames(cdata) %in% samples,]
  }

  rv <- rowVars(mat)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,
                                                     length(rv)))]
  pca <- prcomp(t(mat[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)

  if (!all(intgroup %in% names(cdata))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(cdata[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    cdata[[intgroup]]
  }


  xlab <- paste0("PC", pcx, ": ", round(percentVar[pcx] * 100), "% variance")
  ylab <- paste0("PC", pcy, ": ", round(percentVar[pcy] * 100), "% variance")
  if(is.null(pcz)){
    d <- data.frame(PC1 = pca$x[, pcx],
                    PC2 = pca$x[, pcy],
                    group = group, intgroup.df, name = cdata[,1])


    p <- plot_ly(d,
                 x=~PC1, y=~PC2,
                 color=~group,
                 type='scatter',
                 text=paste0('<b>', d[, 'group'],
                             '</b>\nsample: ', d[, 'name']),
                 mode='markers',
                 hoverinfo='text',
                 marker=list(size=12)) %>%
          layout(xaxis=list(title=xlab,
                            zeroline=FALSE,
                            showline=TRUE,
                            showgrid=FALSE,
                            mirror=TRUE),
                 yaxis=list(title=ylab,
                            zeroline=FALSE,
                            showline=TRUE,
                            showgrid=FALSE,
                            mirror=TRUE))

  } else {
    zlab <- paste0("PC", pcz, ": ", round(percentVar[pcz] * 100), "% variance")

    d <- data.frame(PC1 = pca$x[, pcx],
                    PC2 = pca$x[, pcy],
                    PC3 = pca$x[, pcz],
                    group = group, intgroup.df, name = cdata[,1])

    p <- plot_ly(d,
                 x=~PC1, y=~PC2, z=~PC3,
                 color=~group,
                 type='scatter3d',
                 text=paste0('<b>', d[, 'group'],
                             '</b>\nsample: ', d[, 'name']),
                 mode='markers',
                 hoverinfo='text',
                 marker=list(size=8)) %>%
          layout(scene = list(
                              xaxis=list(title=xlab),
                              yaxis=list(title=ylab),
                              zaxis=list(title=zlab))
                          )

  }

  # add loadings traces
  if(loadings){

    # get loadings for selected PCs
    rot <- pca$rotation

    # columns to use
    cols <- c(pcx, pcy)
    if(!is.null(pcz)) cols <- c(cols, pcz)

    # calculate sum of squares (vector length^2)
    ssq <- rot[, pcx]**2 + rot[, pcy]**2
    if(!is.null(pcz)) ssq <- ssq + rot[, pcz]**2

    # order SS to get top genes
    ssq <- ssq[order(ssq, decreasing=TRUE)]

    top_genes <- names(ssq)[seq_len(loadings_ngenes)]

    # double up a single gene to avoid single-row df effects
    if(length(top_genes) == 1) top_genes <- c(top_genes, top_genes)

    # subset loadings df
    rot_df <- rot[top_genes, ]

    # calculate scale ratio
    all_ratios <- lapply(cols, function(x){
                    diff(range(pca$x[, x]))/diff(range(rot[, x]))
                  })
    scaleratio <- round(mean(unlist(all_ratios)))

    # build text
    tt2 <- lapply(seq_len(nrow(rot_df)), function(x){
             ll <- list(x=rot_df[x, pcx]*scaleratio*1.1,
                        y=rot_df[x, pcy]*scaleratio*1.1,
                        showarrow=FALSE,
                        xanchor='center',
                        yanchor='center',
                        text=rownames(rot_df)[x],
                        font=list(color='red', size=10))
             if(!is.null(pcz)) ll[[ 'z' ]] <- rot_df[x, pcz]*scaleratio*1.1
             ll
           })

    # build line data frame
    line_df <- as.data.frame(rot_df)
    line_df$gene <- rownames(line_df)
    rownames(line_df) <- NULL

    # add origin points
    line_df0 <- line_df
    line_df0[, cols] <- 0
    line_df <- rbind(line_df, line_df0)
    line_df <- unique(line_df)

    if(is.null(pcz)){
      p <- p %>%
              add_trace(x=line_df[, pcx]*scaleratio,
                        y=line_df[, pcy]*scaleratio,
                        type='scatter',
                        mode='lines',
                        line=list(color='black'),
                        split=line_df$gene,
                        showlegend=FALSE,
                        inherit=FALSE) %>%
              layout(annotations=tt2,
                     xaxis=list(zeroline=TRUE),
                     yaxis=list(zeroline=TRUE))
    } else {
      p <- p %>%
              add_trace(x=line_df[, pcx]*scaleratio,
                        y=line_df[, pcy]*scaleratio,
                        z=line_df[, pcz]*scaleratio,
                        type='scatter3d',
                        mode='lines',
                        line=list(color='black'),
                        split=line_df$gene,
                        showlegend=FALSE,
                        inherit=FALSE) %>%
              layout(scene=list(annotations=tt2))
    }

  }

  return(p)
}

#' Make an enrichResult obj from a data frame
#'
#' Most of the parameters are just placeholders and the
#' dataframe must contain the columns 'ID' and 'geneID'
#'
#' @param df data frame with functional enrichment results
#' @param split string, character used to split gene IDs
#' @param keytype type of gene ID
#' @param ontology ontology database being used
#' @param type string, can be 'enrichResult' or 'gseaResult'
#'
#' @return enrichResult object
#'
#' @examples
#' # get enrichResult object
#' data(eres_dex, package='carnation')
#'
#' # extract the results
#' df <- as.data.frame(eres_dex)
#'
#' # convert to a stripped down enrichResult object
#' eres2 <- makeEnrichResult(df)
#'
#' @export
makeEnrichResult <- function(df, split='/',
                             keytype="UNKNOWN",
                             ontology="UNKNOWN",
                             type='enrichResult'){
  if(type == 'enrichResult'){
    if(!all(c('geneID', 'ID') %in% colnames(df))){
      stop('Dataframe must contain the columns "geneID" & "ID"!')
    }

    x <- new('enrichResult',
          result=df,
          keytype=keytype,
          ontology=ontology)
  } else if(type == 'gseaResult'){
    if(!all(c('core_enrichment', 'ID') %in% colnames(df))){
      stop('Dataframe must contain the columns "core_enrichment" & "ID"!')
    }

    x <- new('gseaResult',
          result=df,
          keytype=keytype,
          setType=ontology)
  } else {
    stop('"type" must be "enrichResult" or "gseaResult"!')
  }
  x
}

#' Plot a scatterplot to compare two contrasts
#'
#' @param compare string, what values to plot? can be 'log2FoldChange' or 'P-adj'
#' @param df data frame with log2FoldChange & padj values to plot from 2 contrasts
#' @param label_x string, label for x-axis
#' @param label_y string, label for y-axis
#' @param lim.x x-axis limits
#' @param lim.y y-axis limits
#' @param color.palette character vector of colors to use for significance categories 'Both - same LFC sign',
#'                      'Both - opposite LFC sign', 'None', label_x, label_y
#' @param lab.genes genes to label (default=NULL)
#' @param plot_all string, can be 'yes' or 'no'. if 'yes', points outside axis limits are plotted
#'                 along x/y axis lines (default='no').
#' @param name.col gene name column to merge the 2 results, also used for labeling points
#' @param lines 3-element character vector to plot gridlines in the order (x=0, y=0, x=y),
#'              with 'yes' or 'no' values. E.g. ('yes', 'yes', 'no') will plot dotted lines for
#'              x = 0 & y = 0, but not the x = y diagonal.
#' @param alpha float, marker opacity (default=1).
#' @param size float, marker size (default=4).
#' @param show.grid string, can be 'yes' (default) or 'no'.
#'
#' @return ggplot handle
#'
#' @examples
#' # make mock results df
#' n_genes <- 100
#' res1 <- data.frame(
#'           baseMean = runif(n_genes, 10, 1000),
#'           log2FoldChange = rnorm(n_genes, 0, 2),
#'           lfcSE = runif(n_genes, 0.1, 0.5),
#'           stat = rnorm(n_genes, 0, 3),
#'           pvalue = runif(n_genes, 0, 1),
#'           padj = runif(n_genes, 0, 1),
#'           symbol = paste0("GENE", 1:n_genes),
#'           row.names = paste0("gene", 1:n_genes)
#'         )
#'
#' res2 <- data.frame(
#'           baseMean = runif(n_genes, 10, 1000),
#'           log2FoldChange = rnorm(n_genes, 0, 2),
#'           lfcSE = runif(n_genes, 0.1, 0.5),
#'           stat = rnorm(n_genes, 0, 3),
#'           pvalue = runif(n_genes, 0, 1),
#'           padj = runif(n_genes, 0, 1),
#'           symbol = paste0("GENE", 1:n_genes),
#'           row.names = paste0("gene", 1:n_genes)
#'         )
#'
#' # add geneid column
#' res1 <- cbind(geneid=rownames(res1), res1)
#' res2 <- cbind(geneid=rownames(res2), res2)
#'
#' # make merged df from the two comparisons
#' cols.sub <- c('log2FoldChange', 'padj', 'geneid')
#' df_full <- dplyr::inner_join(
#'   dplyr::select(as.data.frame(res1), all_of(cols.sub)),
#'   dplyr::select(as.data.frame(res2), all_of(cols.sub)),
#'   by = 'geneid',
#'   suffix = c('.x', '.y')
#' )
#'
#' # calculate x & y limits for log2FoldChange
#' xlim <- range(df_full[[ 'log2FoldChange.x' ]])
#' ylim <- range(df_full[[ 'log2FoldChange.y' ]])
#'
#' # get color palette
#' color.palette <- RColorBrewer::brewer.pal(n=5, name='Set2')
#'
#' # add significance column
#' sig.x <- df_full$padj.x < 0.1 & !is.na(df_full$padj.x)
#' sig.y <- df_full$padj.y < 0.1 & !is.na(df_full$padj.y)
#' up.x <- df_full$log2FoldChange.x >= 0
#' up.y <- df_full$log2FoldChange.y >= 0
#' significance <- rep('None', nrow(df_full))
#' significance[ sig.x & sig.y & ((up.x & up.y) | (!up.x & !up.y)) ] <- 'Both - same LFC sign'
#' significance[ sig.x & sig.y & ((up.x & !up.y) | (!up.x & up.y)) ] <- 'Both - opposite LFC sign'
#' significance[ sig.x & !sig.y ] <- 'A vs B'
#' significance[ !sig.x & sig.y ] <- 'B vs A'
#' df_full$significance <- significance
#'
#' # generate scatter plot
#' p <- plotScatter.label(compare = 'log2FoldChange',
#'                        df = df_full,
#'                        label_x = 'A vs B',
#'                        label_y = 'B vs A',
#'                        lim.x = xlim,
#'                        lim.y = ylim,
#'                        color.palette = color.palette)
#'
#' @export
plotScatter.label <- function(compare,
                              df,
                              label_x,
                              label_y,
                              lim.x,
                              lim.y,
                              color.palette,
                              lab.genes=NULL,
                              plot_all='no',
                              name.col='geneid',
                              lines=c('yes', 'yes', 'yes'),
                              alpha=1,
                              size=4,
                              show.grid='yes') {

  names(color.palette) <- c('None', label_x, label_y, 'Both - opposite LFC sign', 'Both - same LFC sign')

  # Determine x and y based on compare
  x <- paste0(compare, '.x')
  y <- paste0(compare, '.y')

  # add '-log10' if padj
  if(compare == 'padj'){
    label_x <- paste('-log10', label_x)
    label_y <- paste('-log10', label_y)
  }

  if (label_x == label_y) {
    label_x <- paste0(label_x, '_x')
    label_y <- paste0(label_y, '_y')
  }

  p <- ggplot(df %>% arrange(.data$significance),
    aes(x = !!sym(x), y = !!sym(y), color = .data$significance, shape = .data$shape, label = !!sym(name.col))) +
    geom_point(alpha=alpha) +
    theme_bw() +
    scale_color_manual(values=color.palette) +
    theme(axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12)) +
    labs(color = 'Significant in:') +
    scale_y_continuous(limits = c(lim.y[1] + 0.02*lim.y[1], lim.y[2] + 0.02*lim.y[2]), expand = c(0.0, 0)) +
    scale_x_continuous(limits = c(lim.x[1] + 0.02*lim.x[1], lim.x[2] + 0.02*lim.x[2]), expand = c(0.0, 0)) +
    xlab(label_x) +
    ylab(label_y)

  if (plot_all == 'yes') {
    p <- p + scale_shape_manual(breaks=c('in','above','below','left','right'),
                                values=c(16, 2, 6, 3, 3),
                                guide='none') +
             scale_size_manual(values=c(size, size+1, size+1, size+1, size+1))
  } else if (plot_all == 'no') {
    p <- p + scale_shape_manual(breaks=c('in'),
                                values=c(16),
                                guide='none') +
             scale_size_manual(values=c(size))
  }
  # Lines
  if (lines[1] == 'yes') {
    p <- p + geom_vline(xintercept=0, color="#333333", linetype="dashed", linewidth=0.5, alpha=0.7)
  }
  if (lines[2] == 'yes') {
    p <- p + geom_hline(yintercept=0, color="#333333", linetype="dashed", linewidth=0.5, alpha=0.7)
  }
   if (lines[3] == 'yes') {
    p <- p + geom_abline(intercept=0, slope=1, color="#333333", linetype="dashed", linewidth=0.5, alpha=0.7)
  }
  # show.grid
  show.grid <- if(show.grid == 'yes') TRUE else FALSE
  # Now apply the theme adjustments
  p <- p + theme(
    panel.grid.major = if(show.grid) element_line() else element_blank(),
    panel.grid.minor = if(show.grid) element_line() else element_blank()
  )

  # Add gene labels if specified
  if (!is.null(lab.genes)) {
    lab.list <- df %>% filter(!!sym(name.col) %in% lab.genes)
    p <- p + geom_point(data=lab.list, color="black", shape=1, size=3)
    p <- p + geom_label_repel(data=lab.list,
      aes(label = !!sym(name.col)),
      fontface="italic",
      show.legend=FALSE,
      max.overlaps=Inf)
    }

  return(p)
}

#' Plot an interactive scatterplot to compare two contrasts
#'
#' @param compare string, what values to plot? can be 'log2FoldChange' or 'P-adj'
#' @param df data frame with log2FoldChange & padj values to plot from 2 contrasts
#' @param label_x string, label for x-axis
#' @param label_y string, label for y-axis
#' @param lim.x x-axis limits
#' @param lim.y y-axis limits
#' @param color.palette character vector of colors to use for significance categories 'Both - same LFC sign',
#'                      'Both - opposite LFC sign', 'None', label_x, label_y
#' @param lab.genes genes to label (default=NULL)
#' @param name.col gene name column to merge the 2 results, also used for labeling points
#' @param lines 3-element character vector to plot gridlines in the order (x=0, y=0, x=y),
#'              with 'yes' or 'no' values. E.g. ('yes', 'yes', 'no') will plot dotted lines for
#'              x = 0 & y = 0, but not the x = y diagonal.
#' @param alpha float, marker opacity (default=1).
#' @param size float, marker size (default=4).
#' @param show.grid string, can be 'yes' (default) or 'no'.

#' @return plotly handle
#'
#' @examples
#' # make mock results df
#' n_genes <- 100
#' res1 <- data.frame(
#'           baseMean = runif(n_genes, 10, 1000),
#'           log2FoldChange = rnorm(n_genes, 0, 2),
#'           lfcSE = runif(n_genes, 0.1, 0.5),
#'           stat = rnorm(n_genes, 0, 3),
#'           pvalue = runif(n_genes, 0, 1),
#'           padj = runif(n_genes, 0, 1),
#'           symbol = paste0("GENE", 1:n_genes),
#'           row.names = paste0("gene", 1:n_genes)
#'         )
#'
#' res2 <- data.frame(
#'           baseMean = runif(n_genes, 10, 1000),
#'           log2FoldChange = rnorm(n_genes, 0, 2),
#'           lfcSE = runif(n_genes, 0.1, 0.5),
#'           stat = rnorm(n_genes, 0, 3),
#'           pvalue = runif(n_genes, 0, 1),
#'           padj = runif(n_genes, 0, 1),
#'           symbol = paste0("GENE", 1:n_genes),
#'           row.names = paste0("gene", 1:n_genes)
#'         )
#'
#' # add geneid column
#' res1 <- cbind(geneid=rownames(res1), res1)
#' res2 <- cbind(geneid=rownames(res2), res2)
#'
#' # make merged df from the two comparisons
#' cols.sub <- c('log2FoldChange', 'padj', 'geneid')
#' df_full <- dplyr::inner_join(
#'   dplyr::select(as.data.frame(res1), all_of(cols.sub)),
#'   dplyr::select(as.data.frame(res2), all_of(cols.sub)),
#'   by = 'geneid',
#'   suffix = c('.x', '.y')
#' )
#'
#' # calculate x & y limits for log2FoldChange
#' xlim <- range(df_full[[ 'log2FoldChange.x' ]])
#' ylim <- range(df_full[[ 'log2FoldChange.y' ]])
#'
#' # get color palette
#' color.palette <- RColorBrewer::brewer.pal(n=5, name='Set2')
#'
#' # add significance column
#' sig.x <- df_full$padj.x < 0.1 & !is.na(df_full$padj.x)
#' sig.y <- df_full$padj.y < 0.1 & !is.na(df_full$padj.y)
#' up.x <- df_full$log2FoldChange.x >= 0
#' up.y <- df_full$log2FoldChange.y >= 0
#' significance <- rep('None', nrow(df_full))
#' significance[ sig.x & sig.y & ((up.x & up.y) | (!up.x & !up.y)) ] <- 'Both - same LFC sign'
#' significance[ sig.x & sig.y & ((up.x & !up.y) | (!up.x & up.y)) ] <- 'Both - opposite LFC sign'
#' significance[ sig.x & !sig.y ] <- 'A vs B'
#' significance[ !sig.x & sig.y ] <- 'B vs A'
#' df_full$significance <- significance
#'
#' # generate scatter plot
#' p <- plotScatter.label_ly(compare = 'log2FoldChange',
#'                           df = df_full,
#'                           label_x = 'A vs B',
#'                           label_y = 'B vs A',
#'                           lim.x = xlim,
#'                           lim.y = ylim,
#'                           color.palette = color.palette)
#'
#' @export
plotScatter.label_ly <- function(compare,
                                 df,
                                 label_x,
                                 label_y,
                                 lim.x,
                                 lim.y,
                                 color.palette,
                                 lab.genes=NULL,
                                 name.col='geneid',
                                 lines=c('yes', 'yes', 'yes'),
                                 alpha=1,
                                 size=4,
                                 show.grid='yes') {

  names(color.palette) <- c('None', label_x, label_y, 'Both - opposite LFC sign', 'Both - same LFC sign')

  # Determine x and y based on compare
  x <- paste0(compare, '.x')
  y <- paste0(compare, '.y')

  label_x <- paste0(compare, ': ', label_x)
  label_y <- paste0(compare, ': ', label_y)

  # add '-log10' if padj
  if(compare == 'padj'){
    label_x <- paste('-log10', label_x)
    label_y <- paste('-log10', label_y)
  }

  if (label_x == label_y) {
    label_x <- paste0(label_x, '_x')
    label_y <- paste0(label_y, '_y')
  }

  show.grid <- if (show.grid == 'yes') TRUE else FALSE

  p <- plot_ly()

  # list of plotting characters
  pch <- c('in'='circle',
           'above'='triangle-up-open',
           'below'='triangle-down-open',
           'left'='triangle-left-open',
           'right'='triangle-right-open')

  # Loop over each level of the factor and create a trace
  for (level in levels(df$significance)) {
    df_i <- df[df$significance == level, ]

    # further subdivide based on shape column
    all_sym <- unique(df_i[['shape']])
    for(sym in all_sym){
      # don't show legend if point is outside limits
      if(sym != 'in'){
        ps <- size + 1

        # only show legend if the *only* points in this level are outside
        if(!all(all_sym == sym)) showlegend <- FALSE
      } else {
        ps <- size
        showlegend <- TRUE
      }
      p <- p %>% add_trace(data = df_i[df_i$shape == sym, ],
                           x = ~get(x),
                           y = ~get(y),
                           type = 'scatter',
                           mode = 'markers',
                           text = ~get(name.col),
                           hoverinfo = 'text',
                           marker = list(size = ps,
                                         opacity = alpha,
                                         symbol = pch[sym],
                                         color = color.palette[level]),
                           showlegend = showlegend,
                           name = level) # This will be the legend entry
    }
  }

  # add labels if any
  if(!is.null(lab.genes)){
    # get data frame of genes to be labeled
    lab.df <- df[df[[name.col]] %in% lab.genes, ]
    if(nrow(lab.df) > 0){
      p <- p %>%
          add_markers(x=lab.df[[x]], y=lab.df[[y]],
                      text=lab.df$name.col,
                      hoverinfo='none',
                      name='Gene scratchpad',
                      marker=list(color='black',
                                  symbol='circle-open',
                                  size=size*2,
                                  line=list(width=2)))
    }
  }

  # Now add the layout
  p <- p %>% layout(xaxis = list(title = label_x,
                                 range = c(lim.x[1] + 0.02*lim.x[1], lim.x[2] + 0.02*lim.x[2]),
                                 zeroline = FALSE,
                                 showgrid = show.grid,
                                 linecolor = 'black',
                                 linewidth = 0.5,
                                 mirror = TRUE),
                    yaxis = list(title = label_y,
                                 range = c(lim.y[1] + 0.02*lim.y[1], lim.y[2] + 0.02*lim.y[2]),
                                 zeroline = FALSE,
                                 showgrid = show.grid,
                                 linecolor = 'black',
                                 linewidth = 0.5,
                                 mirror = TRUE),
                    legend = list(x = 1.05,
                                  y = 0.5,
                                  orientation = 'v',
                                  xanchor = 'left',
                                  yanchor = 'middle',
                                  font = list(size = 14),
                                  title = list(text = 'Significant in:')))

  # Determine the range for the diagonal line
  diag_range_min <- max(lim.x[1], lim.y[1])  # Start from the higher of the two minimums
  diag_range_max <- min(lim.x[2], lim.y[2])  # End at the lower of the two maximums

  # Define shapes for lines including the diagonal line
  shapes_list <- list(
    # Vertical line at x=0
    list(type = "line", x0 = 0, x1 = 0, y0 = 0, y1 = 1, xref = "x", yref = "paper", line = list(color = "grey", dash = "dot")),
    # Horizontal line at y=0
    list(type = "line", x0 = 0, x1 = 1, y0 = 0, y1 = 0, xref = "paper", yref = "y", line = list(color = "grey", dash = "dot")),
    # Diagonal line x=y
    list(type = "line", x0 = diag_range_min, x1 = diag_range_max, y0 = diag_range_min, y1 = diag_range_max, xref = "x", yref = "y", line = list(color = "grey", dash = "dot"))
  )

  # Initialize an empty list for shapes to add
  shapes_to_add <- list()

  if (lines[1] == 'yes') {
    shapes_to_add <- c(shapes_to_add, shapes_list[1])
  }
  if (lines[2] == 'yes') {
    shapes_to_add <- c(shapes_to_add, shapes_list[2])
  }
  if (lines[3] == 'yes') {
    shapes_to_add <- c(shapes_to_add, shapes_list[3])
  }

  # Add all collected shapes at once
  if (length(shapes_to_add) > 0) {
    p <- p %>% layout(shapes = shapes_to_add)
  }

  # Return the Plotly plot object
  return(p)
}

#' Make dummy GeneTonic object
#'
#' @param eres enrichResult object
#'
#' @return GeneTonic object
dummy_genetonic <- function(eres){

  if(!inherits(eres, 'enrichResult')) return(NULL)
  eres2 <- GeneTonic::shake_enrichResult(eres)

  ## add dummy/placeholder columns
  eres2[[ 'z_score' ]] <- 1
  eres2[[ 'aggr_score' ]] <- 0

  # add dummy gene mapping
  all_genes <- unique(unlist(lapply(eres2$gs_genes, function(x) strsplit(x, '\\,')[[1]])))
  anno_df <- data.frame(gene_id=all_genes, gene_name=all_genes, row.names=all_genes)

  return(
    list(l_gs=eres2, anno_df=anno_df)
  )
}

#' Make example carnation object
#'
#' Returns example carnation object used in examples & testing
#'
#' @return reactiveValues object containing carnation object
#'
#' @examplesIf interactive()
#' obj <- make_example_carnation_object()
#'
#' @export
make_example_carnation_object <- function(){
  dds <- makeExampleDESeqDataSet()
  rld <- varianceStabilizingTransformation(dds, blind=TRUE)

  dds <- DESeq(dds)
  results <- results(dds, contrast=c('condition', 'A', 'B'))

  # Create reactive values to simulate app state
  obj <- list(
    dds = list(main = dds),
    rld = list(main = rld),
    res = list(comp1 = results),
    all_dds = dds,
    all_rld = rld,
    dds_mapping = list(comp1 = 'main')
  )

  obj
}

#' Generate upset plot table
#'
#' @param gene.lists list with character vectors of gene names
#' @param comp_split_pattern character used to separate gene set names
#'
#' @return list with upset table elements
#'
#' @examples
#' lst <- list(group1 = c(a = "gene1", b = "gene2", c = "gene3", d = "gene4"),
#'             group2 = c(b = "gene2", d = "gene4", e = "gene5"),
#'             group3 = c(d = "gene4", e = "gene5", f = "gene6"))
#'
#' df <- get_upset_table(lst)
#' str(df)
#'
#' @export
get_upset_table <- function(gene.lists, comp_split_pattern=';'){
  # get matrix of intersections, add set column & save
  df <- fromList.with.names(gene.lists)
  df <- add.set.column(df)

  # get mapping of set names to intersects
  df_unique <- unique(df[, setdiff(colnames(df), 'symbol')])
  rownames(df_unique) <- df_unique$set
  df_unique <- df_unique[, setdiff(colnames(df_unique), 'set')]
  if(nrow(df_unique) == 1){
    set_mapping <- setNames(list(colnames(df_unique)), rownames(df_unique))
  } else {
    set_mapping <- apply(df_unique, 1, function(x){
                     n <- colnames(df_unique)[x == 1]
                     n
        })
  }

  # get columns with degree & comparisons & add to df
  comps <- unlist(
             lapply(set_mapping[df$set],
               function(x) paste(x, collapse=comp_split_pattern)
             )
           )
  degree <- unlist(lapply(set_mapping[df$set], length))

  tbl <- cbind(df, comparisons=comps, degree=degree)

  # get intersection sizes
  inter_counts <- table(df$set)

  # create named vector of choices & add groups + size to names
  # e.g. if 'set1' has 150 genes,
  # then the name is 'set1 (n = 150)'
  inter_choices <- names(inter_counts)
  names(inter_choices) <- paste0(inter_choices, ' (',
                          unlist(
                            lapply(set_mapping,
                              function(x)
                                paste(x, collapse=', ')
                            )
                          ), '; n = ', inter_counts, ')')

  list(tbl=tbl, set_mapping=set_mapping, set_labels=inter_choices)
}

