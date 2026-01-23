# Helper function to create mock DESeq2 data for testing
create_mock_dds <- function(n_genes = 100, n_samples = 6) {
  # round up to even number
  if((n_samples %% 2) != 0) n_samples <- n_samples + 1

  # Create count matrix
  counts <- matrix(
    rpois(n_genes * n_samples, lambda = 100),
    nrow = n_genes,
    ncol = n_samples
  )
  rownames(counts) <- paste0("gene", 1:n_genes)
  colnames(counts) <- paste0("sample", 1:n_samples)

  # Create sample metadata
  coldata <- data.frame(
    condition = factor(rep(c("control", "treatment"), each = n_samples/2)),
    batch = factor(rep(c("A", "B"), times = n_samples/2)),
    samplename = colnames(counts),
    row.names = colnames(counts)
  )

  # Create DESeqDataSet
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = coldata,
    design = ~ condition
  )

  return(dds)
}

# Helper function to create mock results
create_mock_results <- function(n_genes = 100) {
  data.frame(
    baseMean = runif(n_genes, 10, 1000),
    log2FoldChange = rnorm(n_genes, 0, 2),
    lfcSE = runif(n_genes, 0.1, 0.5),
    stat = rnorm(n_genes, 0, 3),
    pvalue = runif(n_genes, 0, 1),
    padj = runif(n_genes, 0, 1),
    symbol = paste0("GENE", 1:n_genes),
    row.names = paste0("gene", 1:n_genes)
  )
}

mock_carnation_obj <- function(){
  dds <- create_mock_dds()
  res <- create_mock_results()
  rld <- DESeq2::varianceStabilizingTransformation(dds, blind=TRUE)

  dds_list <- list(main=dds)
  res_list <- list(comp1=res)
  rld_list <- list(main=rld)
  all_dds <- dds
  all_rld <- rld
  dds_mapping <- list(comp1="main")

  list(
    res=res_list,
    dds=dds_list,
    rld=rld_list,
    all_dds=all_dds,
    all_rld=all_rld,
    dds_mapping=dds_mapping,
    enrich=NULL,
    genetonic=NULL,
    degpatterns=NULL
  )
}

# Helper function to create mock degpatterns data
create_mock_degpatterns <- function(n_genes = 50, n_samples = 6) {
  # Create sample metadata
  coldata <- data.frame(
    condition = factor(rep(c("control", "treatment"), each = n_samples/2)),
    group = factor(rep(c("A", "B", "C"), length.out = n_samples)),
    samplename = paste0("sample", 1:n_samples),
    row.names = paste0("sample", 1:n_samples)
  )

  # Create normalized expression data
  genes <- paste0("gene", 1:n_genes)
  symbols <- paste0("GENE", 1:n_genes)

  # Create clusters (3 clusters)
  clusters <- sample(1:3, n_genes, replace = TRUE)

  # Create expression data for each gene across samples
  expr_data <- lapply(1:n_genes, function(i) {
    # Base expression with cluster-specific pattern
    base_expr <- rnorm(n_samples, mean = 5 + clusters[i], sd = 0.5)

    data.frame(
      genes = genes[i],
      symbol = symbols[i],
      value = base_expr,
      cluster = clusters[i],
      condition = coldata$condition,
      group = coldata$group,
      samplename = coldata$samplename,
      stringsAsFactors = FALSE
    )
  })

  # Combine all genes into one data frame
  normalized <- do.call(rbind, expr_data)

  # Return in degPatterns format
  list(normalized = normalized)
}

# Helper function to create mock enrichment data
create_mock_enrichment <- function(n_terms = 10, n_genes = 20) {
  genes <- paste0("gene", 1:n_genes)

  # make geneID list & get lenth
  geneID <- sapply(1:n_terms, function(i) sample(genes, sample(5:10, 1)))
  len <- unlist(lapply(geneID, length))

  # get background distribution
  bg <- sample(50:150, n_terms, replace = TRUE)
  nbg <- 1000

  # hypergeometric p-values
  pval <- phyper(len, bg, nbg - bg, n_genes, lower.tail=FALSE)
  padj <- p.adjust(pval, method='fdr')

  # Create enrichment result data frame
  df <- data.frame(
    ID = paste0("GO:", sprintf("%07d", 1:n_terms)),
    Description = paste("Biological Process", 1:n_terms),
    GeneRatio = paste0(len, "/", n_genes),
    BgRatio = paste0(bg, "/", nbg),
    pvalue = pval,
    p.adjust = padj,
    qvalue = padj,
    geneID = sapply(geneID, function(x) paste(x, collapse='/')),
    #pvalue = runif(n_terms, 0.0001, 0.05),
    #p.adjust = runif(n_terms, 0.001, 0.1),
    #qvalue = runif(n_terms, 0.001, 0.1),
    #geneID = sapply(1:n_terms, function(i) paste(geneID, collapse = "/")),
    Count = len,
    stringsAsFactors = FALSE
  )

  df
}

# Helper function to create mock genetonic data
create_mock_genetonic <- function(enrich_df, res_df) {
  # Convert geneID from "/" to ","
  l_gs <- enrich_df
  l_gs$gs_id <- l_gs$ID
  l_gs$gs_description <- l_gs$Description
  l_gs$gs_pvalue <- l_gs$pvalue
  l_gs$gs_genes <- gsub("/", ",", l_gs$geneID)
  l_gs$gs_de_count <- l_gs$Count
  l_gs$gs_bg_count <- as.numeric(sapply(strsplit(l_gs$BgRatio, "/"), `[`, 1))
  l_gs$z_score <- rnorm(nrow(l_gs), 0, 1)
  l_gs$aggr_score <- runif(nrow(l_gs), 0, 1)

  # Create annotation data frame
  all_genes <- unique(unlist(strsplit(l_gs$gs_genes, ",")))
  anno_df <- data.frame(
    gene_id = all_genes,
    gene_name = toupper(all_genes),
    row.names = all_genes,
    stringsAsFactors = FALSE
  )

  list(l_gs = l_gs, anno_df = anno_df)
}


