% Generated from roxygen comments by tools/roxy2rd.py
\name{fit_distance_decay}
\alias{fit_distance_decay}
\title{Spatial distance-decay regressions}
\usage{
fit_distance_decay(data, responses, distances, x = "X", y = "Y", correlation = c("exp", "gaus", "spher"), log_distance = TRUE, average_duplicates = TRUE)
}
\arguments{
\item{data}{Data frame with responses, distances and coordinates.}
\item{responses}{Response columns (e.g. factor contributions and
their \code{_logit} relative contributions).}
\item{distances}{Distance columns (e.g. from \code{\link{distance_to_features}()}).}
\item{x,y}{Coordinate columns (projected units).}
\item{correlation}{One of \code{"exp"}, \code{"gaus"}, \code{"spher"}.}
\item{log_distance}{Use \code{log1p(distance)} as the predictor.}
\item{average_duplicates}{See Details.}
}
\description{
Fits, for every response x distance pair, a generalized least-squares
model \code{response ~ log1p(distance)} with a spatial correlation structure
on the sample coordinates (\pkg{nlme}), and returns the slope of the
distance term. A negative slope means the response increases toward
the features.

Duplicate coordinates are averaged first (\code{average_duplicates = TRUE}),
because coincident points make spatial correlation structures singular.
}
\value{
Data frame with \code{response}, \code{distance}, \code{estimate},
\code{std_error}, \code{t_value}, \code{p_value}, \code{n}; failed fits are returned
with \code{NA} and a \code{note}.
}
