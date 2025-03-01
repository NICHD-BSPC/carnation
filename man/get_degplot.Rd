% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{get_degplot}
\alias{get_degplot}
\title{Plot a degPatterns object}
\usage{
get_degplot(
  obj,
  time,
  color = NULL,
  cluster_column = "cluster",
  cluster_to_show,
  x_order,
  points = TRUE,
  boxes = TRUE,
  smooth = "smooth",
  lines = TRUE,
  facet = TRUE,
  prefix_title = "Cluster ",
  genes_to_label = NULL
)
}
\arguments{
\item{obj}{degPatterns object}

\item{time}{metadata variable to plot on x-axis}

\item{color}{variable to color plot}

\item{cluster_column}{column to use for grouping genes}

\item{cluster_to_show}{which clusters to show in plot}

\item{x_order}{order of x-axis values}

\item{points}{boolean, show samples on plot? Default: TRUE}

\item{boxes}{boolean, show boxes on plot? Default: TRUE}

\item{smooth}{what type of trendline to use? can be 'smooth' (default) or 'line'.}

\item{lines}{show lines joining samples? Default: TRUE}

\item{facet}{boolean, should plot be faceted? Default: TRUE}

\item{prefix_title}{string, prefix for facet titles}

\item{genes_to_label}{genes to label on plot}
}
\value{
ggplot handle
}
\description{
This function plots a degPatterns object using
degPlotCluster
}
