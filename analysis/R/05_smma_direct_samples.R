# SMMA applied directly to raw samples, bypassing PMF (Section 3.6,
# Supplementary Section S11, Table S10).
source("analysis/R/00_config.R")

samples <- read.csv(file.path(cfg$data_dir, "pmf_input_samples.csv"), check.names = FALSE)
srclib <- read_source_library(cfg$ex("source_library.csv"))
kow     <- read_kow(cfg$ex("kow.csv"))
catmap  <- read.csv(file.path(cfg$data_dir, "source_category_mapping.csv"), check.names = FALSE)
groups  <- analyte_groupings(names(samples))
sites   <- read.csv(file.path(cfg$out_dir, "site_classification.csv"))   # from script 02

# Label library profiles with the nine main-text categories, so per-sample
# category summaries are in that scheme (the fits do not depend on labels).
srclib$category <- catmap$category_maintext[match(srclib$source, catmap$source)]

# EPA16 composition only, Kow weights, no S/N (S/N is a PMF diagnostic)
epa16 <- intersect(names(samples), groups$EPA16)    # column order of the input file
bt <- smma_samples(samples, srclib, sample_col = "pit_code", B = 300, seed = cfg$seed,
                   analytes = epa16, min_shared = 8, min_analytes = 8,
                   kow = kow, kow_args = cfg$kow_args)
out(bt$per_source, "smma_direct_samples_per_source.csv")
out(bt$per_category, "smma_direct_samples_per_category.csv")

# Per-sample weights summed into the nine main-text categories, then median
# and mean across the samples of each site class.
by_site_class <- function(w) {
  ps <- stats::aggregate(weight ~ factor + category, data = w, FUN = sum)
  names(ps)[1] <- "sample"
  ps$site_class <- sites$site_class[match(ps$sample, sites$sample)]
  s <- do.call(rbind, lapply(split(ps, list(ps$site_class, ps$category), drop = TRUE),
                             function(d) data.frame(site_class = d$site_class[1],
                                                    category = d$category[1],
                                                    n = length(unique(d$sample)),
                                                    median_pct = 100 * stats::median(d$weight),
                                                    mean_pct = 100 * mean(d$weight))))
  list(per_sample = ps,
       summary = s[order(match(s$site_class, c("Urban", "Mixed", "Remote")), -s$median_pct), ])
}

# Table S10: unweighted single NNLS fit per sample (point estimate; does not
# depend on the Kow weights).
pt <- smma_samples(samples, srclib, sample_col = "pit_code", B = 0,
                   analytes = epa16, min_shared = 8, min_analytes = 8)
s10 <- by_site_class(pt)
out(s10$per_sample, "smma_direct_samples_point_by_sample.csv")
out(s10$summary, "smma_direct_samples_table_s10.csv")
print(s10$summary, row.names = FALSE)

# Bootstrapped best-guess mixture per sample (Kow-weighted; Figure S21).
# Medians and means are over the samples in which the category was retained.
mix <- smma_best_mixture(bt, min_sel_rate = cfg$min_sel_rate, min_weight = cfg$min_weight)
bs <- by_site_class(mix$sources)
out(bs$per_sample, "smma_direct_samples_best_mixture_by_sample.csv")
out(bs$summary, "smma_direct_samples_bootstrap_by_site_class.csv")
