% Generated from roxygen comments by tools/roxy2rd.py
\name{read_sn}
\alias{read_sn}
\title{Read factor-specific signal-to-noise (S/N) values}
\usage{
read_sn(x, sheet = "SN", factor_col = "factor")
}
\arguments{
\item{x}{File path or data frame with one row per factor and one
column per analyte (as in the \code{SN} sheet written for EPA PMF output).}
\item{sheet}{Excel sheet (default \code{"SN"}).}
\item{factor_col}{Name of the factor-identifier column.}
}
\description{
Read factor-specific signal-to-noise (S/N) values
}
\value{
A data frame.
}
