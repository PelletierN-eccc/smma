% Generated from roxygen comments by tools/roxy2rd.py
\name{read_factor_contributions}
\alias{read_factor_contributions}
\title{Read receptor-model factor contributions (sample scores)}
\usage{
read_factor_contributions(x, sheet = "Contributions", sample_col = "sample")
}
\arguments{
\item{x}{A file path (CSV or Excel) or a data frame with one row per
sample and one column per factor.}
\item{sheet}{Excel sheet (default \code{"Contributions"}).}
\item{sample_col}{Name of the sample-identifier column.}
}
\description{
Read receptor-model factor contributions (sample scores)
}
\value{
A data frame with \code{sample_col} and numeric factor columns.
}
