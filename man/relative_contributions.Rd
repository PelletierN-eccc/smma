% Generated from roxygen comments by tools/roxy2rd.py
\name{relative_contributions}
\alias{relative_contributions}
\title{Relative factor contributions and their logits}
\usage{
relative_contributions(contributions, factors, eps = 0.001)
}
\arguments{
\item{contributions}{Data frame of factor contributions.}
\item{factors}{Factor columns.}
\item{eps}{Bounds used before the logit transform.}
}
\description{
Relative factor contributions and their logits
}
\value{
\code{contributions} with \code{<factor>_pct} and \code{<factor>_logit}
columns added (negative contributions set to 0 first).
}
