% Generated from roxygen comments by tools/roxy2rd.py
\name{smma_best_mixture}
\alias{smma_best_mixture}
\title{Best-guess source mixture from bootstrap results}
\usage{
smma_best_mixture(boot, min_sel_rate = 0.20, min_weight = 0.01, top_n = NULL, stat = c("median", "mean"))
}
\arguments{
\item{boot}{An \code{smma_bootstrap} object.}
\item{min_sel_rate,min_weight}{Retention thresholds.}
\item{top_n}{Optional maximum number of sources per factor.}
\item{stat}{Statistic used as the weight: \code{"median"} or \code{"mean"}.}
}
\description{
Keeps, within each factor, the sources whose selection rate is at
least \code{min_sel_rate} or whose median (or mean) weight is at least
\code{min_weight} (optionally only the \code{top_n} largest), and renormalizes
their weights to sum to 1. The top-ranked source is always kept.
}
\value{
List with \code{sources} (factor, source, category, weight and the
bootstrap summary columns) and \code{categories} (weight summed by category).
}
\examples{
ex <- smma_example()
bt <- smma_bootstrap(smma_prepare(ex$profiles, ex$library), B = 50, seed = 1)
smma_best_mixture(bt)$categories
}
