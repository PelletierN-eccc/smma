% Generated from roxygen comments by tools/roxy2rd.py
\name{update_library}
\alias{update_library}
\title{Replace library values for selected sources}
\usage{
update_library(source_library, values, source_col = "source", analytes = NULL, analyte_pattern = "^(pah_|tph_|biom_)")
}
\arguments{
\item{source_library}{Source library.}
\item{values}{Data frame with the source identifier column and the
analyte columns to overwrite (other columns are ignored).}
\item{source_col}{Source identifier column (in both tables).}
\item{analytes}{Analyte columns to take from \code{values} (default: those
matching \code{analyte_pattern}).}
\item{analyte_pattern}{Pattern used when \code{analytes} is \code{NULL}.}
}
\description{
Overwrites analyte values of the sources listed in \code{values} (one row
per source, one column per analyte to replace). Typical use: build
library variants with bounding values for analytes that literature
profiles do not report, then refit to test sensitivity.
}
\value{
The modified library.
}
