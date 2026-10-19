% Generated from roxygen comments by tools/roxy2rd.py
\name{kow_weights}
\alias{kow_weights}
\title{Per-analyte partitioning (Kow) weights}
\usage{
kow_weights(analytes, kow = NULL, input_scale = c("log10", "linear"), transform = NULL, matrix = c("sediment", "water"))
}
\arguments{
\item{analytes}{Character vector of analyte names.}
\item{kow}{Named numeric vector (see \code{\link{read_kow}()}); \code{NULL} gives all 1.}
\item{input_scale}{Scale of the values in \code{kow}: \code{"log10"} if they are
already log10 Kow (the usual tabulated form), \code{"linear"} if they are
Kow. Weights are proportional to log10 Kow in both cases.}
\item{transform}{Optional function applied to the log10 Kow values
before scaling (e.g. \code{log10} reproduces analyses in which tabulated
log Kow values were log-transformed a second time).}
\item{matrix}{\code{"sediment"} (weight increases with Kow) or \code{"water"}
(reversed, so more soluble compounds weigh more).}
}
\description{
Builds one weight per analyte from octanol-water partition
coefficients so that particle-affinitive compounds weigh more in a
sediment matrix (or less, for water). Weights are scaled to a mean of 1
over the analytes that have a usable value; analytes without one get 1.
}
\value{
Named numeric vector of weights.
}
\examples{
kow <- read_kow(system.file("extdata", "kow.csv", package = "smma"))
head(kow_weights(names(kow), kow))
}
