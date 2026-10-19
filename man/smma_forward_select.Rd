% Generated from roxygen comments by tools/roxy2rd.py
\name{smma_forward_select}
\alias{smma_forward_select}
\title{Forward selection of a small number of sources}
\usage{
smma_forward_select(input, max_k = 3, threshold = 0.10, min_analytes = 12)
}
\arguments{
\item{input}{An \code{smma_input} object from \code{\link{smma_prepare}()}.}
\item{max_k}{Maximum number of sources.}
\item{threshold}{Relative RMSE at which selection stops.}
\item{min_analytes}{Minimum number of usable analytes per factor.}
}
\description{
Greedy selection of up to \code{max_k} library sources per factor, adding at
each step the source that most reduces the relative RMSE, and stopping
once the relative RMSE falls below \code{threshold}.
}
\value{
Data frame with \code{factor}, \code{step}, \code{source}, \code{category},
\code{weight} (final weights of the selected set) and \code{rel_rmse} (after
each step).
}
