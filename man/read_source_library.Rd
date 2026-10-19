% Generated from roxygen comments by tools/roxy2rd.py
\name{read_source_library}
\alias{read_source_library}
\title{Read a source-profile library}
\usage{
read_source_library(x, source_col = "source", category_col = "category")
}
\arguments{
\item{x}{File path or data frame with one row per source profile, an
identifier column, a category column and one column per analyte.}
\item{source_col,category_col}{Identifier and category column names.}
}
\description{
Read a source-profile library
}
\value{
A data frame (source identifiers are made unique if needed).
}
