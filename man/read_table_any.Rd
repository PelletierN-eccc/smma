% Generated from roxygen comments by tools/roxy2rd.py
\name{read_table_any}
\alias{read_table_any}
\title{Read a table from CSV or Excel}
\usage{
read_table_any(path, sheet = NULL)
}
\arguments{
\item{path}{File path.}
\item{sheet}{Sheet name or number for Excel files (ignored for CSV).}
}
\description{
Small convenience reader used by the other \code{read_*()} functions. CSV
files are read with \code{utils::read.csv()} (a UTF-8 byte-order mark is
handled); \code{.xlsx}/\code{.xls} files require the \pkg{readxl} package.
}
\value{
A data frame.
}
