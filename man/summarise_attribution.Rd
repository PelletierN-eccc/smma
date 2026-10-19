% Generated from roxygen comments by tools/roxy2rd.py
\name{summarise_attribution}
\alias{summarise_attribution}
\title{Summarise attributed concentrations by group}
\usage{
summarise_attribution(attribution, groups = NULL, group_col = "group", level = "category", B = 2000, conf = 0.95, seed = NULL)
}
\arguments{
\item{attribution}{Output of \code{\link{attribute_sources}()}.}
\item{groups}{Data frame mapping \code{sample} to a grouping column, or
\code{NULL} to summarise all samples together.}
\item{group_col}{Name of the grouping column in \code{groups}.}
\item{level}{Column of \code{attribution} holding categories (or sources).}
\item{B}{Bootstrap resamples for the confidence interval (0 to skip).}
\item{conf}{Confidence level.}
\item{seed}{Optional seed.}
}
\description{
Computes, for each group (e.g. urban/mixed/remote site class) x
category x analyte grouping, the number of samples, the mean
attributed concentration with a percentile-bootstrap confidence
interval (resampling samples), the median, and the share of the
group's mean burden.
}
\value{
Data frame with \code{group}, \code{level}, \code{grouping}, \code{n}, \code{mean},
\code{lo}, \code{hi}, \code{median}, \code{share_pct}.
}
