% Generated from roxygen comments by tools/roxy2rd.py
\name{analyte_groupings}
\alias{analyte_groupings}
\title{Define analyte groupings for source attribution}
\usage{
analyte_groupings(analytes, scheme = ests_analyte_scheme())
}
\arguments{
\item{analytes}{Character vector of available analyte column names
(e.g. \code{names(profiles)}).}
\item{scheme}{Analyte scheme, see \code{\link{ests_analyte_scheme}()}.}
}
\description{
Returns the four groupings used in the accompanying manuscript
(EPA16, Parent, Alkylated, Total) for a given set of analyte columns.
Any named list of character vectors can be used instead in
\code{\link{attribute_sources}()}.
}
\value{
A named list of character vectors.
}
\examples{
prof <- read_factor_profiles(system.file("extdata", "bc_factor_profiles.csv",
package = "smma"))
str(analyte_groupings(names(prof)))
}
