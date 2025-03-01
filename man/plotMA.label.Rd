% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{plotMA.label}
\alias{plotMA.label}
\title{Create a labeled MA plot}
\usage{
plotMA.label(
  res,
  fdr.thres = 0.01,
  fc.thres = 0,
  fc.lim = NULL,
  lab.genes = NULL,
  tolower.cols = c("SYMBOL", "ALIAS")
)
}
\arguments{
\item{res}{data.frame with DE analysis results. Must contain
"padj" & "log2FoldChange" columns}

\item{fdr.thres}{False discovery rate (FDR) threshold}

\item{fc.thres}{log2FoldChange threshold}

\item{fc.lim}{y-axis limits}

\item{lab.genes}{genes to label on MA plot}

\item{tolower.cols}{column names that will be converted to
lower case}
}
\value{
ggplot handle
}
\description{
This function creates an MA plot from a data.frame
containing DE analysis results.
}
