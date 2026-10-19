# Manuscript tables (Table 2; Tables S4, S6, S10, S12, S13, S14) as CSV and
# formatted .xlsx, laid out as in the article. Needs the output of scripts
# 01, 02 and 04-06. Writes to <out_dir>/tables/.
source("analysis/R/00_config.R")

tab_dir <- file.path(cfg$out_dir, "tables")
dir.create(tab_dir, showWarnings = FALSE)
rd <- function(f) utils::read.csv(file.path(cfg$out_dir, f), check.names = FALSE)
dash <- function(x) gsub("Medium-heavy", "Medium–heavy", gsub(" - ", " – ", x))
f1 <- function(x) formatC(x, format = "f", digits = 1, big.mark = ",")
site_lev <- c("Urban", "Mixed", "Remote")
basis_lev <- c("EPA16", "Parent", "Total", "Alkylated")

write_table <- function(x, name, title, note = NULL) {
  utils::write.csv(x, file.path(tab_dir, paste0(name, ".csv")), row.names = FALSE)
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "table")
  openxlsx::writeData(wb, "table", title, startRow = 1)
  openxlsx::addStyle(wb, "table", openxlsx::createStyle(textDecoration = "bold"), rows = 1, cols = 1)
  openxlsx::writeData(wb, "table", x, startRow = 3,
                      headerStyle = openxlsx::createStyle(textDecoration = "bold",
                                                          border = "TopBottom", wrapText = TRUE))
  if (!is.null(note)) openxlsx::writeData(wb, "table", note, startRow = nrow(x) + 5)
  openxlsx::setColWidths(wb, "table", cols = seq_along(x),
                         widths = pmin(60, pmax(10, nchar(names(x)) + 2,
                                                vapply(x, function(v) max(nchar(as.character(v)), 0), 1))))
  openxlsx::freezePane(wb, "table", firstActiveRow = 4)
  openxlsx::saveWorkbook(wb, file.path(tab_dir, paste0(name, ".xlsx")), overwrite = TRUE)
}

# --- Table 2: mean (95% CI); share of the site-class mean burden --------------
t2 <- rd("attribution_mean_ci_by_site_class_9cat.csv")
t2 <- t2[!grepl("^Light petroleum", t2$category), ]
t2$cell <- sprintf("%s (%s–%s); %s%%", f1(t2$mean), f1(t2$lo), f1(t2$hi), f1(t2$share_pct))
w <- stats::reshape(t2[, c("group", "category", "n", "grouping", "cell")], direction = "wide",
                    idvar = c("group", "category", "n"), timevar = "grouping")
names(w) <- sub("^cell\\.", "", names(w))
epa <- t2[t2$grouping == "EPA16", ]
w$order <- epa$mean[match(paste(w$group, w$category), paste(epa$group, epa$category))]
w <- w[order(match(w$group, site_lev), -w$order), c("group", "category", "n", basis_lev)]
w$category <- dash(w$category)
names(w)[1:2] <- c("Site category", "Source category")
write_table(w, "Table_2_attribution_by_site_class",
            "Table 2. Mean attributed concentration (ng g-1 dw, 95% CI) and percentage of the site-category mean burden, by analyte basis and site category.",
            "Values are mean attributed concentration (95% CI); percentage of the site-category mean burden. Light petroleum received no attribution and is omitted.")

# --- Table S4: Tukey outliers --------------------------------------------------
s4 <- rd("tukey_outliers_table_s4.csv")
s4 <- data.frame(Sample = s4$sample,
                 `Sum16 EPA PAC (ng g-1)` = formatC(s4$sum16, format = "f", digits = 0, big.mark = ","),
                 `Largest attributed category (%)` = sprintf("%s (%s%%)", dash(s4$first), f1(s4$first_pct)),
                 `Second (%)` = sprintf("%s (%s%%)", dash(s4$second), f1(s4$second_pct)),
                 `All petroleum (%)` = f1(s4$petroleum_pct), check.names = FALSE)
write_table(s4, "Table_S4_tukey_outliers",
            "Table S4. Samples exceeding the Tukey outlier threshold for log10-transformed sum of 16 EPA PACs, with their largest attributed source categories.")

# --- Table S6: share versus % fines -----------------------------------------------
s6 <- rd("share_vs_fines_table_s6.csv")
s6 <- data.frame(`Source category` = dash(s6$category),
                 `Pearson r vs % fines` = sprintf("%+.3f", s6$r),
                 `p value` = formatC(s6$p_value, format = "g", digits = 2),
                 `Significant (p < 0.05)` = ifelse(s6$p_value < 0.05, "Yes", "No"),
                 check.names = FALSE)
write_table(s6, "Table_S6_share_vs_fines",
            sprintf("Table S6. Correlation between source-category percentage contribution (EPA16 basis) and percent fine sediment (n = %d).",
                    rd("share_vs_fines_table_s6.csv")$n[1]))

# --- Table S10: SMMA applied directly to samples (point estimate) ----------------
s10 <- rd("smma_direct_samples_table_s10.csv")
s10 <- data.frame(`Site category` = s10$site_class, `Source category` = dash(s10$category),
                  n = s10$n, `Median %` = f1(s10$median_pct), `Mean %` = f1(s10$mean_pct),
                  check.names = FALSE)
write_table(s10, "Table_S10_direct_smma",
            "Table S10. Median and mean percentage contribution of each source category to EPA16 burden, by site category, from a single unweighted NNLS fit of each raw sample (SMMA applied directly to samples).")

# --- Table S12: biomarker sensitivity -------------------------------------------
s12 <- rd("biomarker_sensitivity_petrogenic_share.csv")
s12 <- s12[order(match(s12$magnitude, c("baseline", "low", "central", "high")),
                 s12$capped_profiles_excluded, s12$factor), ]
s12 <- data.frame(`Biomarker magnitude bound` = tools::toTitleCase(s12$magnitude),
                  `Capped-ratio profiles excluded` = ifelse(s12$capped_profiles_excluded, "Yes", "No"),
                  Factor = s12$factor, `Median petrogenic share (%)` = f1(100 * s12$median),
                  `95% CI (%)` = sprintf("%s–%s", f1(100 * s12$p025), f1(100 * s12$p975)),
                  `Sources in library` = s12$n_sources, check.names = FALSE)
write_table(s12, "Table_S12_biomarker_sensitivity",
            "Table S12. Median bootstrap petrogenic share and 95% CI for Factor 1 and Factor 2 across four biomarker-magnitude bounds, with and without the three capped-ratio profiles.")

# --- Table S13: median per-sample shares, 17 categories -------------------------
s13 <- rd("attribution_median_share_by_site_class_17cat.csv")
w <- stats::reshape(s13[, c("group", "category", "n", "grouping", "median_pct")], direction = "wide",
                    idvar = c("group", "category", "n"), timevar = "grouping")
names(w) <- sub("^median_pct\\.", "", names(w))
w <- w[rowSums(w[, basis_lev]) > 0, ]
w <- w[order(match(w$group, site_lev), -w$EPA16), c("group", "category", "n", basis_lev)]
w[basis_lev] <- lapply(w[basis_lev], f1)
w$category <- dash(w$category)
names(w)[1:2] <- c("Site category", "Source category")
write_table(w, "Table_S13_median_shares_17cat",
            "Table S13. Median per-sample percentage contribution of each source category, 17-category detailed scheme, by analyte basis and site category.")

# --- Table S14: category-level SMMA weights per factor -------------------------
s14 <- rd("smma_best_mixture_17cat.csv")
s14 <- s14[order(s14$factor, -s14$weight), ]
s14 <- data.frame(Factor = s14$factor, `Source category` = dash(s14$category),
                  `SMMA weight` = paste0(f1(100 * s14$weight), "%"), check.names = FALSE)
write_table(s14, "Table_S14_factor_weights_17cat",
            "Table S14. Category-level SMMA best-guess weight per PMF factor, 17-category detailed scheme (companion to main-text Figure 3).")
cat("Tables written to", tab_dir, "\n")
