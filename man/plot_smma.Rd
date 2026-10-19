% Generated from roxygen comments by tools/roxy2rd.py
\name{plot_smma}
\alias{plot_smma}
\title{Plot SMMA category weights per factor}
\usage{
plot_smma(boot, min_weight = 0.005)
}
\arguments{
\item{boot}{An \code{smma_bootstrap} object, or the \code{categories} element of
\code{\link{smma_best_mixture}()}.}
\item{min_weight}{Hide categories whose median weight is below this value.}
}
\description{
Plot SMMA category weights per factor
}
\value{
A ggplot object (requires \pkg{ggplot2}).
}
