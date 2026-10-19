% Generated from roxygen comments by tools/roxy2rd.py
\name{default_population_tiers}
\alias{default_population_tiers}
\title{Default population tiers for site classification}
\usage{
default_population_tiers()
}
\description{
Each settlement is assigned a tier by population, and each tier has an
"urban" radius (site within the settlement's built-up/traffic
footprint) and a wider "mixed" radius (partial influence). Tiers
without an urban radius can only produce "Mixed" classifications.
These are the values used for coastal British Columbia in the
accompanying study; adapt them to your region.
}
\value{
Data frame with \code{tier}, \code{min_pop}, \code{urban_km}, \code{mixed_km}.
}
