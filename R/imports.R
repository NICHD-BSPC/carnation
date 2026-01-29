#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom colorspace rainbow_hcl qualitative_hcl
#' @importFrom ComplexUpset upset upset_modify_themes intersection_size upset_query
#' @importFrom dendextend as.ggdend get_branches_heights
#' @importFrom DESeq2 normalizationFactors sizeFactors design counts DESeqDataSetFromMatrix
#' estimateSizeFactors plotPCA varianceStabilizingTransformation makeExampleDESeqDataSet DESeq results
#' @importFrom dplyr mutate relocate select filter any_of "%>%" rename inner_join case_when all_of
#' arrange everything
#' @importFrom DT renderDT DTOutput datatable formatStyle formatSignif dataTableProxy selectRows selectCells
#' @importFrom enrichplot cnetplot
#' @importFrom GeneTonic cluster_markov distill_enrichment gs_horizon
#' gs_summary_overview gs_summary_overview_pair enrichment_map gs_dendro
#' gs_alluvial gs_fuzzyclustering shake_enrichResult get_aggrscores shake_gsenrichResult
#' @import ggplot2
#' @importFrom ggrepel geom_label_repel
#' @importFrom graphics par
#' @importFrom grDevices dev.off pdf colorRampPalette
#' @importFrom heatmaply heatmaply is.plotly
#' @importFrom htmltools withTags tagAppendChild tagAppendChildren tags tagList
#' @importFrom igraph V "V<-"
#' @importFrom MatrixGenerics rowVars
#' @importFrom methods new
#' @importFrom plotly plotlyOutput renderPlotly layout plot_ly add_trace add_markers toWebGL save_image ggplotly
#' @importFrom RColorBrewer brewer.pal
#' @importFrom reticulate py_install use_virtualenv
#' @importFrom rintrojs introjsUI introBox introjs
#' @importFrom scales alpha
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @importFrom shinyBS bsCollapse bsCollapsePanel updateCollapse
#' @importFrom shinycssloaders withSpinner
#' @importFrom shinymanager secure_app secure_server create_db check_credentials
#' @importFrom shinythemes shinytheme
#' @importFrom shinyWidgets dropdownButton tooltipOptions
#' @importFrom sortable bucket_list add_rank_list
#' @importFrom stats as.formula prcomp setNames median
#' @importFrom SummarizedExperiment "colData<-" colData assay
#' @importFrom tools file_path_sans_ext
#' @importFrom utils read.table packageName
#' @importFrom viridisLite viridis
#' @importFrom visNetwork renderVisNetwork visNetworkOutput visIgraph visOptions visNodes
#' @importFrom yaml read_yaml write_yaml
NULL
