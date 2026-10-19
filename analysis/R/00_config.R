# =========================================================================
# Configuration for reproducing Pelletier et al. (submitted, Environmental
# Pollution): "Source Apportionment of PACs and Petroleum Biomarkers in
# British Columbia Intertidal Sediments: An Urban-to-Remote Baseline for
# Spill Response".
#
# Run the scripts from the repository root, in order:
#   source("analysis/R/01_smma_factors.R")
#   source("analysis/R/02_source_attribution.R")
#   source("analysis/R/03_spatial_regression.R")      # needs DFO layers, see README
#   source("analysis/R/04_biomarker_sensitivity.R")   # slow (8 x B runs)
#   source("analysis/R/05_smma_direct_samples.R")
# =========================================================================

library(smma)   # remotes::install_git("https://gccode.ssc-spc.gc.ca/PelletierN-ECCC/smma")

cfg <- list(
  data_dir   = "analysis/data",
  ex         = function(f) system.file("extdata", f, package = "smma", mustWork = TRUE),

  # SMMA settings used in the manuscript
  B          = getOption("smma.B", 10000),  # bootstrap draws per factor
  seed       = 42,
  min_shared = 12,
  min_sel_rate = 0.20,  # best-guess mixture retention rules
  min_weight   = 0.01,

  # Kow weighting. kow.csv tabulates log10 Kow, and weights are
  # proportional to log10 Kow, scaled to a mean of 1 (Methods 2.6).
  # kow_double_log = TRUE reproduces the weighting of the original
  # submission draft, in which these values were log10-transformed a
  # second time (weights proportional to log10(log10 Kow)).
  kow_double_log = getOption("smma.kow_double_log", FALSE),

  # Attribution summaries
  ci_B = 2000, ci_seed = 42
)
# Outputs of the two weightings are kept apart.
cfg$out_dir <- if (cfg$kow_double_log) "analysis/output_double_log" else "analysis/output"
cfg$kow_args <- if (cfg$kow_double_log) {
  list(input_scale = "log10", transform = log10, matrix = "sediment")
} else {
  list(input_scale = "log10", matrix = "sediment")
}
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
out <- function(x, f) utils::write.csv(x, file.path(cfg$out_dir, f), row.names = FALSE)
