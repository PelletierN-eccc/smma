# SMMA applied directly to raw samples, bypassing PMF (Section 3.6,
# Supplementary Section S11).
source("analysis/R/00_config.R")

samples <- read.csv(file.path(cfg$data_dir, "pmf_input_samples.csv"), check.names = FALSE)
srclib <- read_source_library(cfg$ex("source_library.csv"))
kow     <- read_kow(cfg$ex("kow.csv"))
catmap  <- read.csv(file.path(cfg$data_dir, "source_category_mapping.csv"), check.names = FALSE)
groups  <- analyte_groupings(names(samples))

# EPA16 composition only, Kow weights, no S/N (S/N is a PMF diagnostic)
bt <- smma_samples(samples, srclib, sample_col = "pit_code", B = 300, seed = cfg$seed,
                   analytes = groups$EPA16, min_shared = 10, min_analytes = 10,
                   kow = kow, kow_args = cfg$kow_args)
bt$per_source$category_maintext <- catmap$category_maintext[match(bt$per_source$source, catmap$source)]
out(bt$per_source, "smma_direct_samples_per_source.csv")
out(bt$per_category, "smma_direct_samples_per_library_category.csv")
