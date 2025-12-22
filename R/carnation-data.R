#' A `DESeqResults` object testing the effect of dexamethasone
#' on smooth muscle cells
#'
#' @format A `DESeqResults` object, generated in the `DESeq2` framework
#'
#' @details This `DESeqResults` object on the data from the `airway` package
#' has been created comparing dexamethasone treated vs untreated samples,
#' accounting for the different cell lines included.
#'
#' Details on how this object has been created are included in the
#' `create_carnation_data.R` script, included in the `scripts` folder of the
#' `Carnation` package.
#'
#' @references Himes BE, Jiang X, Wagner P, Hu R, Wang Q, Klanderman B, Whitaker
#' RM, Duan Q, Lasky-Su J, Nikolos C, Jester W, Johnson M, Panettieri R Jr,
#' Tantisira KG, Weiss ST, Lu Q. “RNA-Seq Transcriptome Profiling Identifies
#' CRISPLD2 as a Glucocorticoid Responsive Gene that Modulates Cytokine Function
#' in Airway Smooth Muscle Cells.” PLoS One. 2014 Jun 13;9(6):e99625. PMID:
#' 24926665. GEO: GSE52778
#'
#' @name res_dex
#' @docType data
NULL

#' A `DESeqResults` object testing the difference between
#' two cell lines of smooth muscle cells
#'
#' @format A `DESeqResults` object, generated in the `DESeq2` framework
#'
#' @details This `DESeqResults` object on the data from the `airway` package
#' has been created comparing two smooth muscle cell lines,
#' accounting for the effect of dexamethasone treatment.
#'
#' Details on how this object has been created are included in the
#' `create_carnation_data.R` script, included in the `scripts` folder of the
#' `Carnation` package.
#'
#' @references Himes BE, Jiang X, Wagner P, Hu R, Wang Q, Klanderman B, Whitaker
#' RM, Duan Q, Lasky-Su J, Nikolos C, Jester W, Johnson M, Panettieri R Jr,
#' Tantisira KG, Weiss ST, Lu Q. “RNA-Seq Transcriptome Profiling Identifies
#' CRISPLD2 as a Glucocorticoid Responsive Gene that Modulates Cytokine Function
#' in Airway Smooth Muscle Cells.” PLoS One. 2014 Jun 13;9(6):e99625. PMID:
#' 24926665. GEO: GSE52778
#'
#' @name res_cell
#' @docType data
NULL

#' An `enrichResult` object for differentially expressed genes in the
#' dexamethasone treatment comparison.
#'
#' @format An `enrichResult` object, generated with the `enrichGO` function
#' from the `clusterProfiler` package.
#'
#' @details This `enrichResult` object was created to test for functional
#' enrichment using the GO Biological Process (BP) ontology on the top
#' 100 differentially expressed genes from the dexamethasone treatment
#' comparison.
#'
#' Details on how this object has been created are included in the
#' `create_carnation_data.R` script, included in the `scripts` folder of the
#' `Carnation` package.
#'
#' @references Himes BE, Jiang X, Wagner P, Hu R, Wang Q, Klanderman B, Whitaker
#' RM, Duan Q, Lasky-Su J, Nikolos C, Jester W, Johnson M, Panettieri R Jr,
#' Tantisira KG, Weiss ST, Lu Q. “RNA-Seq Transcriptome Profiling Identifies
#' CRISPLD2 as a Glucocorticoid Responsive Gene that Modulates Cytokine Function
#' in Airway Smooth Muscle Cells.” PLoS One. 2014 Jun 13;9(6):e99625. PMID:
#' 24926665. GEO: GSE52778
#'
#' @name enrich_bp_dex
#' @docType data
NULL

#' An `enrichResult` object for differentially expressed genes in the
#' cell line comparison.
#'
#' @format An `enrichResult` object, generated with the `enrichGO` function
#' from the `clusterProfiler` package.
#'
#' @details This `enrichResult` object was created to test for functional
#' enrichment using the GO Biological Process (BP) ontology on the top
#' 100 differentially expressed genes from the cell line
#' comparison.
#'
#' Details on how this object has been created are included in the
#' `create_carnation_data.R` script, included in the `scripts` folder of the
#' `Carnation` package.
#'
#' @references Himes BE, Jiang X, Wagner P, Hu R, Wang Q, Klanderman B, Whitaker
#' RM, Duan Q, Lasky-Su J, Nikolos C, Jester W, Johnson M, Panettieri R Jr,
#' Tantisira KG, Weiss ST, Lu Q. “RNA-Seq Transcriptome Profiling Identifies
#' CRISPLD2 as a Glucocorticoid Responsive Gene that Modulates Cytokine Function
#' in Airway Smooth Muscle Cells.” PLoS One. 2014 Jun 13;9(6):e99625. PMID:
#' 24926665. GEO: GSE52778
#'
#' @name enrich_bp_cell
#' @docType data
NULL

#' A `degPatterns` object for differentially expressed genes in the
#' dexamethasone treatment comparison.
#'
#' @format A `degPatterns` object, generated with the `degPatterns` function
#' from the `DEGreport` package.
#'
#' @details This `degPatterns` object was created to test for groups of
#' coexpressed genes in the top 100 differentially expressed genes from the
#' dexamethasone treatment comparison.
#'
#' Details on how this object has been created are included in the
#' `create_carnation_data.R` script, included in the `scripts` folder of the
#' `Carnation` package.
#'
#' @references Himes BE, Jiang X, Wagner P, Hu R, Wang Q, Klanderman B, Whitaker
#' RM, Duan Q, Lasky-Su J, Nikolos C, Jester W, Johnson M, Panettieri R Jr,
#' Tantisira KG, Weiss ST, Lu Q. “RNA-Seq Transcriptome Profiling Identifies
#' CRISPLD2 as a Glucocorticoid Responsive Gene that Modulates Cytokine Function
#' in Airway Smooth Muscle Cells.” PLoS One. 2014 Jun 13;9(6):e99625. PMID:
#' 24926665. GEO: GSE52778
#'
#' @name degpatterns_dex
#' @docType data
NULL


