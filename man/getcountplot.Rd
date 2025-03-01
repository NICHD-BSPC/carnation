% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{getcountplot}
\alias{getcountplot}
\title{Create gene plot}
\usage{
getcountplot(
  df,
  intgroup = "group",
  factor.levels,
  title = NULL,
  ylab = "Normalized counts",
  color = "gene",
  nrow = 2,
  ymin = NULL,
  ymax = NULL,
  log = TRUE,
  freey = FALSE,
  trendline = "smooth",
  facet = NULL,
  legend = TRUE,
  boxes = TRUE,
  rotate_x_labels = 30
)
}
\arguments{
\item{df}{data.frame with gene counts}

\item{intgroup}{metadata variable to plot on x-axis}

\item{factor.levels}{levels of intgroup to show on x-axis}

\item{title}{title of plot}

\item{ylab}{y-axis label}

\item{color}{metadata variable to color by}

\item{nrow}{number of rows to plot if faceting}

\item{ymin}{y-axis lower limit}

\item{ymax}{y-axis upper limit}

\item{log}{should y-axis be log10-transformed?}

\item{freey}{should y-axes of faceted plots have independent scales?}

\item{trendline}{type of trendline to draw}

\item{facet}{metadata variable to facet by}

\item{legend}{show legend?}

\item{boxes}{show boxes?}

\item{rotate_x_labels}{angle to rotate x-axis labels (default=30)}
}
\value{
ggplot handle
}
\description{
This function creates the gene plot.
}
