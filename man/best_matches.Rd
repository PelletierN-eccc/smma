% Generated from roxygen comments by tools/roxy2rd.py
\name{best_matches}
\alias{best_matches}
\title{Best-matching library sources per factor and metric}
\usage{
best_matches(distances, n = 1)
}
\arguments{
\item{distances}{Output of \code{\link{profile_distances}()}.}
\item{n}{Number of matches to return per factor and metric.}
}
\description{
Best-matching library sources per factor and metric
}
\value{
Long data frame: \code{factor}, \code{metric}, \code{rank}, \code{source},
\code{category}, \code{score}, \code{n_shared}.
}
