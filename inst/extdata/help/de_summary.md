#### Summary Table
------------------

Overview table showing differential expression statistics for all comparisons in your analysis.

**What it shows:**
- Number of up and downregulated genes per comparison
- Total genes analyzed and filtering statistics
- Experimental design information

**When to use it:**
- To get a quick overview of your differential expression results
- To compare the number of DE genes across different contrasts
- To verify analysis parameters and filtering steps

**How to interpret:**
- `comparison`: Unique identifier for each differential expression comparison
- `up`, `down`: Number of genes significantly up or downregulated
- `total.genes`: Total number of annotated genes included in the analysis
- `total.nonzero`: Number of genes with non-zero counts in at least one sample
- `outliers`: Genes flagged as count outliers by DESeq2 (high variability between replicates)
- `low.counts`: Genes filtered out due to low expression across all samples
- `design`: Linear model formula used for the DESeq2 analysis

**Note:** Changes to the FDR and LFC thresholds affect all DE analysis visualizations.

