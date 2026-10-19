% Generated from roxygen comments by tools/roxy2rd.py
\name{analyte_scheme}
\alias{analyte_scheme}
\title{Build a custom analyte scheme}
\usage{
analyte_scheme(pac_pattern, biomarker_pattern = NULL, tph_pattern = NULL, alkylated_pattern = NULL, exclude_pattern = NULL, epa16 = character(), ratio = character(), homologues = list(), labels = NULL)
}
\arguments{
\item{pac_pattern,biomarker_pattern,tph_pattern}{Regular expressions
identifying PAC, biomarker and TPH columns.}
\item{alkylated_pattern}{Regular expression identifying alkylated
homologues among the PAC columns.}
\item{exclude_pattern}{Regular expression for summary columns to drop
(use \code{NULL} for none).}
\item{epa16}{Character vector of EPA 16 priority PAH column names.}
\item{ratio}{Named character vector; see \code{\link{ests_analyte_scheme}()}.}
\item{homologues}{Named list of homologue series (parent first).}
\item{labels}{Optional named character vector of display labels.}
}
\description{
Build a custom analyte scheme
}
\value{
A list of class \code{smma_scheme}.
}
