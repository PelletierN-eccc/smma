% Generated from roxygen comments by tools/roxy2rd.py
\name{smma_group_share}
\alias{smma_group_share}
\title{Bootstrap share of a group of sources}
\usage{
smma_group_share(boot, groups, categories = NULL)
}
\arguments{
\item{boot}{An \code{smma_bootstrap} object created with \code{keep_draws = TRUE}.}
\item{groups}{Named character vector mapping source identifiers to groups,
or a function of the category returning the group.}
\item{categories}{Data frame \code{source}, \code{category} (e.g. \code{input$categories}),
required when \code{groups} is a function.}
}
\description{
Sums, within every bootstrap draw, the weights of the sources assigned
to each group (e.g. petrogenic versus other) and summarises the summed
weight with the same fit-weighted statistics as \code{\link{smma_bootstrap}()}.
}
\value{
Data frame \code{factor}, \code{group}, \code{median}, \code{p025}, \code{p975}, \code{mean}.
}
