% Generated from roxygen comments by tools/roxy2rd.py
\name{sn_weights}
\alias{sn_weights}
\title{Factor-specific signal-to-noise weights}
\usage{
sn_weights(sn, factor_id, analytes, factor_col = "factor")
}
\arguments{
\item{sn}{Data frame of S/N values (see \code{\link{read_sn}()}) or \code{NULL}.}
\item{factor_id}{Factor identifier.}
\item{analytes}{Analyte names.}
\item{factor_col}{Factor column in \code{sn}.}
}
\description{
Extracts the S/N row of one factor as a weight vector over \code{analytes}.
Analytes absent from the table get weight 1.
}
\value{
Named numeric vector.
}
