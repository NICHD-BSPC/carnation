library(airway)
library(DESeq2)
library(org.Hs.eg.db)

# load airway data
data('airway')

# extract counts and metadata
mat <- assay(airway)
cdata <- colData(airway)

# filter some genes with low counts in many samples
keep <- rowSums(mat >= 10) >= 6
mat <- mat[keep, ]

# get symbol annotations
anno_df <- mapIds(org.Hs.eg.db,
               column='SYMBOL',
               keys=rownames(mat),
               keytype='ENSEMBL')

# analyze with DESeq2
dds <- DESeqDataSetFromMatrix(mat,
                              colData=cdata,
                              design=~cell + dex)
dds <- DESeq(dds)

# extract comparisons of interest
res_list <- list(
              cell = results(dds, contrast = c("cell", "N61311","N052611")),
              dex = results(dds, contrast = c("dex", "trt", "untrt"))
            )

for(comp in names(res_list)){
  res <- res_list[[ comp ]]

  # add gene column from rownames
  res$gene <- rownames(res)

  # add symbol column from annotations
  res$symbol <- anno_df[rownames(res)]

  # get DE genes with FDR < 0.1
  de_genes <- rownames(res)[res$padj < 0.1 & !is.na(res$padj)]

  # functional enrichment of top genes using GO BP
  eres <- clusterProfiler::enrichGO(
              gene = de_genes[1:100],
              keyType = 'ENSEMBL',
              OrgDb=org.Hs.eg.db,
              ont='BP'
          )

  # save res & enrich objects
  save(res, file = paste0("data/res_", comp, ".RData"), compress="xz")
  save(eres, file = paste0("data/enrich_bp_", comp, ".RData"), compress="xz")
}

# get degPatterns object
ma.i <- mat[rownames(mat) %in% de_genes[1:100],]

# remove any genes with 0 variance
ma.i <- ma.i[rowVars(ma.i) != 0, ]

# run pattern analysis
p <- DEGreport::degPatterns(
        ma.i,
        cdata,
        time='dex',
        plot=FALSE
        )

# save degPatterns object
save(p, file = paste0("data/degpatterns_dex.RData"), compress="xz")

