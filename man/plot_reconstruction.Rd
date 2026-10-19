% Generated from roxygen comments by tools/roxy2rd.py
\name{plot_reconstruction}
\alias{plot_reconstruction}
\title{Plot observed versus reconstructed factor composition}
\usage{
plot_reconstruction(input, mixture, factor_id)
}
\arguments{
\item{input}{An \code{smma_input} object.}
\item{mixture}{Data frame with \code{factor}, \code{source}, \code{weight}
(e.g. \code{smma_best_mixture()$sources}).}
\item{factor_id}{Factor to plot.}
}
\description{
Plot observed versus reconstructed factor composition
}
\value{
A ggplot object.
}
