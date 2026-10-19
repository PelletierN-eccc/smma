# Old (double log10 Kow, original submission draft) versus corrected (log10
# Kow) weighting for every weight-dependent number in the article. Run
# scripts 01-02 and 04-06 with both settings of cfg$kow_double_log first.
# Writes analysis/output/kow_weighting_comparison.csv.
old_dir <- "analysis/output_double_log"
new_dir <- "analysis/output"
rd <- function(dir, f) utils::read.csv(file.path(dir, f), check.names = FALSE)
both <- function(f) list(old = rd(old_dir, f), new = rd(new_dir, f))
rows <- list()
add <- function(location, item, old, new, digits = 1) {
  rows[[length(rows) + 1]] <<- data.frame(location = location, item = item,
                                          old = round(old, digits), new = round(new, digits))
}

# --- Factor mixtures, 17 categories (Figure 3, Table S14, Sections 3.1.1-3.1.4)
m <- both("smma_best_mixture_17cat.csv")
x <- merge(m$old, m$new, by = c("factor", "category"), all = TRUE, suffixes = c("_old", "_new"))
x[is.na(x)] <- 0
x <- x[pmax(x$weight_old, x$weight_new) >= 0.0005, ]
x <- x[order(x$factor, -x$weight_new), ]
for (i in seq_len(nrow(x))) add("Fig. 3 / Table S14 / 3.1", paste(x$factor[i], x$category[i], "(%)"),
                                100 * x$weight_old[i], 100 * x$weight_new[i])

# --- Table 2: mean (CI) and share, nine categories
t2 <- both("attribution_mean_ci_by_site_class_9cat.csv")
x <- merge(t2$old, t2$new, by = c("group", "category", "grouping"), suffixes = c("_old", "_new"))
x <- x[order(match(x$group, c("Urban", "Mixed", "Remote")), x$category,
             match(x$grouping, c("EPA16", "Parent", "Total", "Alkylated"))), ]
for (i in seq_len(nrow(x))) {
  lab <- paste(x$group[i], x$category[i], x$grouping[i])
  add("Table 2", paste(lab, "mean"), x$mean_old[i], x$mean_new[i])
  add("Table 2", paste(lab, "CI low"), x$lo_old[i], x$lo_new[i])
  add("Table 2", paste(lab, "CI high"), x$hi_old[i], x$hi_new[i])
  add("Table 2", paste(lab, "share (%)"), x$share_pct_old[i], x$share_pct_new[i])
}

# --- Section 3.3 statements derived from Table 2 (EPA16 basis unless stated)
derive <- function(t) {
  e <- t[t$grouping == "EPA16", ]
  s <- function(g, k) sum(e$share_pct[e$group == g & e$category %in% k])
  tot <- function(g, basis = "EPA16") sum(t$mean[t$group == g & t$grouping == basis])
  bio <- "Biomass combustion"; veh <- "Vehicle combustion (diesel/gasoline exhaust)"
  petro <- grep("petroleum|Other oils", unique(t$category), value = TRUE)
  pet <- function(g, basis) sum(t$share_pct[t$group == g & t$grouping == basis & t$category %in% petro])
  cl <- c("Urban", "Mixed", "Remote")
  c(urban_over_remote_EPA16 = tot("Urban") / tot("Remote"),
    stats::setNames(sapply(cl, function(g) s(g, c(bio, veh))), paste0("bio_plus_vehicle_", cl)),
    stats::setNames(sapply(cl, function(g) s(g, bio)), paste0("biomass_", cl)),
    stats::setNames(sapply(cl, function(g) s(g, veh)), paste0("vehicle_", cl)),
    stats::setNames(sapply(cl, function(g) s(g, "Street dust")), paste0("street_dust_", cl)),
    stats::setNames(sapply(cl, function(g) s(g, "Industrial/coal combustion")), paste0("industrial_", cl)),
    stats::setNames(sapply(cl, function(g) pet(g, "EPA16")), paste0("petroleum_EPA16_", cl)),
    stats::setNames(sapply(cl, function(g) pet(g, "Alkylated")), paste0("petroleum_Alkylated_", cl)),
    stats::setNames(sapply(cl, function(g) pet(g, "Alkylated") / pet(g, "EPA16")),
                    paste0("petroleum_alk_over_EPA16_", cl)))
}
d_old <- derive(t2$old); d_new <- derive(t2$new)
for (k in names(d_old)) add("3.3 / Abstract / Highlights", k, d_old[[k]], d_new[[k]], 2)

# --- Table S13: median per-sample shares (nine and 17 categories)
for (tag in c("9cat", "17cat")) {
  s <- both(sprintf("attribution_median_share_by_site_class_%s.csv", tag))
  x <- merge(s$old, s$new, by = c("group", "category", "grouping"), suffixes = c("_old", "_new"))
  x <- x[order(match(x$group, c("Urban", "Mixed", "Remote")), -x$median_pct_new), ]
  for (i in seq_len(nrow(x))) add(sprintf("Table S13 / 3.3 (%s)", tag),
                                  paste(x$group[i], x$category[i], x$grouping[i], "median (%)"),
                                  x$median_pct_old[i], x$median_pct_new[i])
}

# --- Table S12 / 3.1.6: petrogenic share across biomarker bounds
s <- both("biomarker_sensitivity_petrogenic_share.csv")
x <- merge(s$old, s$new, by = c("magnitude", "capped_profiles_excluded", "factor"),
           suffixes = c("_old", "_new"))
for (i in seq_len(nrow(x))) {
  lab <- paste(x$factor[i], x$magnitude[i], if (x$capped_profiles_excluded[i]) "cap3 excluded" else "all")
  add("Table S12 / 3.1.6", paste(lab, "median (%)"), 100 * x$median_old[i], 100 * x$median_new[i])
  add("Table S12 / 3.1.6", paste(lab, "CI low (%)"), 100 * x$p025_old[i], 100 * x$p025_new[i])
  add("Table S12 / 3.1.6", paste(lab, "CI high (%)"), 100 * x$p975_old[i], 100 * x$p975_new[i])
}
for (f in unique(x$factor)) {
  y <- x[x$factor == f, ]
  add("3.1.6", paste(f, "median range min (%)"), 100 * min(y$median_old), 100 * min(y$median_new))
  add("3.1.6", paste(f, "median range max (%)"), 100 * max(y$median_old), 100 * max(y$median_new))
  cap <- function(v) {
    a <- y[!y$capped_profiles_excluded, ]; b <- y[y$capped_profiles_excluded, ]
    b <- b[match(a$magnitude, b$magnitude), ]
    100 * max(abs(a[[v]] - b[[v]]))
  }
  add("3.1.6", paste(f, "max change from excluding cap3 (pp)"), cap("median_old"), cap("median_new"), 3)
}

# --- SMMA applied directly to samples: Table S10 (unweighted point estimate)
# and the bootstrapped, Kow-weighted best-guess mixtures (Figure S21)
direct <- c(`Table S10 / 3.6` = "smma_direct_samples_table_s10.csv",
            `Fig. S21 bootstrap / 3.6` = "smma_direct_samples_bootstrap_by_site_class.csv")
for (loc in names(direct)) {
  s <- both(direct[[loc]])
  x <- merge(s$old, s$new, by = c("site_class", "category"), all = TRUE, suffixes = c("_old", "_new"))
  x <- x[order(match(x$site_class, c("Urban", "Mixed", "Remote")), -x$median_pct_new), ]
  for (i in seq_len(nrow(x))) {
    add(loc, paste(x$site_class[i], x$category[i], "median (%)"), x$median_pct_old[i], x$median_pct_new[i])
    add(loc, paste(x$site_class[i], x$category[i], "mean (%)"), x$mean_pct_old[i], x$mean_pct_new[i])
  }
}

# --- Sections 3.4-3.5, Tables S4 and S6
s <- both("share_vs_sum16_spearman.csv")
x <- merge(s$old, s$new, by = "category", suffixes = c("_old", "_new"))
for (i in seq_len(nrow(x))) {
  add("3.4", paste(x$category[i], "Spearman rho vs sum16"), x$rho_old[i], x$rho_new[i], 3)
  add("3.4", paste(x$category[i], "p value"), x$p_value_old[i], x$p_value_new[i], 6)
}
s <- both("tukey_outliers_table_s4.csv")
x <- merge(s$old, s$new, by = "sample", suffixes = c("_old", "_new"))
for (i in seq_len(nrow(x))) {
  add("Table S4 / 3.4", paste(x$sample[i], "largest:", x$first_new[i],
                              if (x$first_old[i] != x$first_new[i]) "(CHANGED)" else ""),
      x$first_pct_old[i], x$first_pct_new[i])
  add("Table S4 / 3.4", paste(x$sample[i], "second:", x$second_new[i],
                              if (x$second_old[i] != x$second_new[i]) "(CHANGED)" else ""),
      x$second_pct_old[i], x$second_pct_new[i])
}
s <- both("tukey_outliers_elevation.csv")
x <- merge(s$old, s$new, by = "category", suffixes = c("_old", "_new"))
for (i in seq_len(nrow(x))) add("3.4", paste(x$category[i], "outlier / rest median ratio"),
                                x$ratio_old[i], x$ratio_new[i], 0)
s <- both("top5pct_largest_source.csv")
for (k in sort(unique(c(s$old$largest, s$new$largest))))
  add("3.4", paste("top 5% samples led by", k, "(count)"), sum(s$old$largest == k), sum(s$new$largest == k), 0)
s <- both("share_vs_fines_table_s6.csv")
x <- merge(s$old, s$new, by = "category", suffixes = c("_old", "_new"))
for (i in seq_len(nrow(x))) {
  add("Table S6 / 3.5", paste(x$category[i], "Pearson r vs fines"), x$r_old[i], x$r_new[i], 3)
  add("Table S6 / 3.5", paste(x$category[i], "p value"), x$p_value_old[i], x$p_value_new[i], 4)
}

res <- do.call(rbind, rows)
res$diff <- res$new - res$old
utils::write.csv(res, file.path(new_dir, "kow_weighting_comparison.csv"), row.names = FALSE)
cat(nrow(res), "rows written to", file.path(new_dir, "kow_weighting_comparison.csv"), "\n")
