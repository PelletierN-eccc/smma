% Generated from roxygen comments by tools/roxy2rd.py
\name{smma_example}
\alias{smma_example}
\title{Example data: British Columbia intertidal sediments}
\usage{
smma_example()
}
\description{
Loads the example files shipped in \code{inst/extdata}: the four-factor
EPA PMF 5.0 solution for 211 intertidal sediment samples from coastal
British Columbia (profiles, contributions and S/N), the 371-profile
PAC/biomarker source library, and log10 Kow values. Field data are
published in the ECCC Data Catalogue ("ESTS Shoreline surveys data").
}
\value{
A list with \code{profiles}, \code{contributions}, \code{sn}, \code{library}, \code{kow}.
}
\examples{
ex <- smma_example()
names(ex)
}
