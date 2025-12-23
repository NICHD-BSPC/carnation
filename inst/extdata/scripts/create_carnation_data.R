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

# cell comparison
res_cell <- res_list[[ 'cell' ]]

# add gene column from rownames
res_cell$gene <- rownames(res_cell)

# add symbol column from annotations
res_cell$symbol <- anno_df[rownames(res_cell)]

# get DE genes with FDR < 0.01 & |log2FoldChange| >= 1
de_genes <- rownames(res_cell)[res_cell$padj < 0.01 & abs(res_cell$log2FoldChange) >= 1 & !is.na(res_cell$padj)]

# functional enrichment of top genes using GO BP
eres_cell <- clusterProfiler::enrichGO(
               gene = de_genes,
               keyType = 'ENSEMBL',
               OrgDb=org.Hs.eg.db,
               ont='BP',
               pvalueCutoff=1,
               qvalueCutoff=1
             )

# rm geneSets slot to reduce obj size
eres_cell@geneSets <- list()

# save res_cell & enrich objects
save(res_cell, file = paste0("data/res_cell.RData"), compress="xz")
save(eres_cell, file = paste0("data/eres_cell.RData"), compress="xz")

# dex comparison
res_dex <- res_list[[ 'dex' ]]

# add gene column from rownames
res_dex$gene <- rownames(res_dex)

# add symbol column from annotations
res_dex$symbol <- anno_df[rownames(res_dex)]

# get DE genes with FDR < 0.01 & |log2FoldChange| >= 1
de_genes <- rownames(res_dex)[res_dex$padj < 0.01 & abs(res_dex$log2FoldChange) >= 1 & !is.na(res_dex$padj)]

# functional enrichment of top genes using GO BP
eres_dex <- clusterProfiler::enrichGO(
              gene = de_genes[1:100],
              keyType = 'ENSEMBL',
              OrgDb=org.Hs.eg.db,
              ont='BP',
              pvalueCutoff=1,
              qvalueCutoff=1
            )

# rm geneSets slot to reduce obj size
eres_dex@geneSets <- list()

# save res_dex & enrich objects
save(res_dex, file = paste0("data/res_dex.RData"), compress="xz")
save(eres_dex, file = paste0("data/eres_dex.RData"), compress="xz")

# get degPatterns object
ma.i <- mat[rownames(mat) %in% de_genes,]

# remove any genes with 0 variance
ma.i <- ma.i[rowVars(ma.i) != 0, ]

# run pattern analysis
degpatterns_dex <- DEGreport::degPatterns(
                     ma.i,
                     cdata,
                     time='cell',
                     col='dex',
                     reduce=TRUE,
                     plot=FALSE
                   )

# save degPatterns object
save(degpatterns_dex, file = paste0("data/degpatterns_dex.RData"), compress="xz")

