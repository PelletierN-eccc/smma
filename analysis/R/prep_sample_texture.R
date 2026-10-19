# Builds analysis/data/sample_texture.csv (% fines < 63 um and TOC per
# sample) from the ESTS Shoreline surveys data (ECCC Data Catalogue, SQLite
# release). Only needed to regenerate that file; requires DBI and RSQLite.
#   Rscript analysis/R/prep_sample_texture.R path/to/PreSCAT_ESTS_1.0.sqlite
#
# % fines: differential particle-size distribution summed over bins with an
# upper diameter <= 63 um for each sub-sample (est_code), averaged over the
# sub-samples of a shoreline segment. Grain size is reported per segment, so
# every sample of a segment gets the same value.
# TOC: variable "total_organic_carbon", matched on pit_code (mean of
# replicates).
args <- commandArgs(trailingOnly = TRUE)
db <- if (length(args)) args[1] else stop("Give the path to the SQLite file.")
con <- DBI::dbConnect(RSQLite::SQLite(), db, flags = RSQLite::SQLITE_RO)
on.exit(DBI::dbDisconnect(con))
gs <- DBI::dbGetQuery(con, "select segment_code, est_code, particle_diameter_max_um as dmax,
                            distribution from grain_size_distribution")
toc <- DBI::dbGetQuery(con, "select pit_code, value from samples
                             where variable = 'total_organic_carbon'")

samples <- read.csv("analysis/data/pmf_input_samples.csv", check.names = FALSE)$pit_code

gs$distribution <- as.numeric(gs$distribution)
gs$dmax <- suppressWarnings(as.numeric(gs$dmax))
per_est <- aggregate(distribution ~ segment_code + est_code, data = gs[gs$dmax <= 63, ], FUN = sum)
fines <- aggregate(distribution ~ segment_code, data = per_est, FUN = mean)
seg <- sub("_[A-Z]+_[0-9]+$", "", samples)            # BURRA_2018_BEL_UIZ_1 -> BURRA_2018_BEL
toc$value <- as.numeric(toc$value)
toc <- aggregate(value ~ pit_code, data = toc, FUN = mean)

out <- data.frame(sample = samples,
                  fines_pct = fines$distribution[match(seg, fines$segment_code)],
                  toc_pct = toc$value[match(samples, toc$pit_code)])
utils::write.csv(out, "analysis/data/sample_texture.csv", row.names = FALSE)
cat(sprintf("%d samples: %d with fines, %d with TOC, %d with both\n", nrow(out),
            sum(!is.na(out$fines_pct)), sum(!is.na(out$toc_pct)),
            sum(!is.na(out$fines_pct) & !is.na(out$toc_pct))))
