% Generated from roxygen comments by tools/roxy2rd.py
\name{read_kow}
\alias{read_kow}
\title{Read octanol-water partition coefficients}
\usage{
read_kow(x, analyte_col = "analyte", value_col = "kow")
}
\arguments{
\item{x}{File path or data frame.}
\item{analyte_col,value_col}{Column names for the long layout.}
}
\description{
Accepts either a one-row "wide" table (one column per analyte) or a
two-column "long" table (\code{analyte}, \code{value}).
}
\value{
A named numeric vector (names are analyte codes).
}
