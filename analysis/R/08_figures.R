# Journal figures (Figure 3; Figures S4-S15, S24-S26). Needs the output of
# scripts 01-02. Writes 600-dpi LZW TIFF (17 cm wide) and PNG previews to
# analysis/figures/ (analysis/figures_double_log/ for the old weighting).
source("analysis/R/00_config.R")
library(ggplot2)

fig_dir <- if (cfg$kow_double_log) "analysis/figures_double_log" else "analysis/figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
save_fig <- function(p, name, height_cm, width_cm = 17) {
  ggsave(file.path(fig_dir, paste0(name, ".tiff")), p, width = width_cm, height = height_cm,
         units = "cm", dpi = 600, device = "tiff", compression = "lzw", bg = "white")
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = width_cm, height = height_cm,
         units = "cm", dpi = 150, bg = "white")
}

res      <- readRDS(file.path(cfg$out_dir, "smma_factors.rds"))
profiles <- read_factor_profiles(cfg$ex("bc_factor_profiles.csv"))
sn       <- read_sn(cfg$ex("bc_factor_sn.csv"))
kow      <- read_kow(cfg$ex("kow.csv"))
catmap   <- read.csv(file.path(cfg$data_dir, "source_category_mapping.csv"), check.names = FALSE)

# --- One colour per nine-category family, same order in every figure --------
# (validated categorical palette, adjacent-pair CVD and normal-vision checks pass)
family_levels <- c("Biomass combustion", "Vehicle combustion (diesel/gasoline exhaust)",
                   "Street dust", "Industrial/coal combustion",
                   "Medium-heavy petroleum - fresh", "Medium-heavy petroleum - weathered",
                   "Other oils (lubricating)", "Light petroleum (diesel)")
family_cols <- stats::setNames(c("#eb6834", "#2a78d6", "#1baf7a", "#eda100",
                                 "#e87ba4", "#008300", "#6250d6", "#e34948"), family_levels)
family_of <- function(cat9) sub("^Light petroleum \\(diesel\\).*", "Light petroleum (diesel)", cat9)
dash <- function(x) gsub("Medium-heavy", "Medium–heavy", gsub(" - ", " – ", x))  # en dashes
fam_labels <- stats::setNames(dash(family_levels), family_levels)

# 17-category order: by family, then the order below within a family
detail <- unique(catmap[, c("category_detailed", "category_maintext")])
detail$family <- family_of(detail$category_maintext)
detail <- detail[order(match(detail$family, family_levels), detail$category_detailed), ]
detail_levels <- detail$category_detailed
family_of_detail <- stats::setNames(detail$family, detail$category_detailed)

theme_pub <- theme_bw(base_size = 8, base_family = "sans") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3),
        strip.background = element_rect(fill = "grey95", colour = "grey70"),
        strip.text = element_text(face = "bold"), legend.key.size = unit(3, "mm"),
        axis.ticks = element_line(colour = "grey50", linewidth = 0.3))
fill_family <- scale_fill_manual(values = family_cols, labels = fam_labels, drop = TRUE,
                                 name = NULL)
pct_lab <- function(x) ifelse(x >= 0.1, sprintf("%.1f%%", x), "<0.1%")

# --- Figure 3: best-guess mixture per factor, 17 categories ----------------
m17 <- read.csv(file.path(cfg$out_dir, "smma_best_mixture_17cat.csv"))
m17 <- m17[m17$weight >= 0.0005, ]
m17$pct <- 100 * m17$weight
m17$family <- factor(family_of_detail[m17$category], family_levels)
m17$category <- factor(m17$category, rev(detail_levels))
p3 <- ggplot(m17, aes(x = .data$pct, y = .data$category, fill = .data$family)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = pct_lab(.data$pct)), hjust = -0.15, size = 2.3, colour = "grey15") +
  facet_wrap(~ factor, ncol = 2, scales = "free_y") +
  scale_x_continuous(limits = c(0, 105), breaks = seq(0, 100, 25), expand = c(0, 0),
                     labels = function(x) paste0(x, "%")) +
  scale_y_discrete(labels = dash) + fill_family +
  labs(x = "SMMA best-guess weight (bootstrap median, renormalized)", y = NULL) +
  theme_pub + theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 3))
save_fig(p3, "Figure_3_SMMA_best_mixture_17cat", 13)

# --- Figure S4: diagnostic ratios -----------------------------------------------
dr <- read.csv(file.path(cfg$out_dir, "factor_diagnostic_ratios.csv"))
dr <- stats::reshape(dr, direction = "long", varying = names(dr)[-1], v.names = "value",
                     timevar = "ratio", times = names(dr)[-1], idvar = "factor")
ratio_lab <- c(Ant_AntPhe = "Ant/(Ant+Phe)", Fla_FlaPyr = "Fla/(Fla+Pyr)",
               BaA_BaAChr = "BaA/(BaA+Chr)", IcdP_IcdPBghiP = "IcdP/(IcdP+BghiP)")
dr$ratio <- factor(ratio_lab[dr$ratio], rev(ratio_lab))
heat_text <- function(v) ifelse(v > 0.6, "white", "grey10")
pS4 <- ggplot(dr, aes(x = .data$factor, y = .data$ratio, fill = .data$value)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.2f", .data$value), colour = heat_text(.data$value)), size = 2.8) +
  scale_colour_identity() +
  scale_fill_gradient(low = "#cde2fb", high = "#104281", limits = c(0, 1), name = "Ratio") +
  labs(x = NULL, y = NULL) + theme_pub + theme(panel.grid = element_blank())
save_fig(pS4, "Figure_S4_diagnostic_ratios", 6.5, 12)

# --- Figures S5-S7: alkylation indices ------------------------------------------
ai <- read.csv(file.path(cfg$out_dir, "factor_alkylation_indices.csv"))
ai$series <- factor(ai$series, rev(names(default_analyte_scheme()$homologues)))
alk_fig <- function(col, legend, name, fmt, trans = "identity", limits = NULL) {
  d <- ai[, c("factor", "series", col)]; names(d)[3] <- "value"
  lab <- ifelse(is.finite(d$value), sprintf(fmt, d$value), "n/a")
  if (trans == "log10") d$value <- pmax(d$value, 0.1)    # zero ratios drawn at the scale floor
  dark <- is.finite(d$value) & d$value > (if (trans == "log10") 5 else 0.55)
  p <- ggplot(d, aes(x = .data$factor, y = .data$series, fill = .data$value)) +
    geom_tile(colour = "white", linewidth = 0.8) +
    geom_text(aes(label = lab), colour = ifelse(dark, "white", "grey10"), size = 2.6) +
    scale_fill_gradient(low = "#cde2fb", high = "#104281", trans = trans, limits = limits,
                        na.value = "grey85", name = legend) +
    labs(x = NULL, y = NULL) + theme_pub + theme(panel.grid = element_blank())
  save_fig(p, name, 8, 13)
}
alk_fig("parent_fraction", "Parent\nfraction", "Figure_S5_parent_fraction", "%.2f", limits = c(0, 1))
alk_fig("alkyl_parent_ratio", "Alkylated/\nparent", "Figure_S6_alkyl_parent_ratio", "%.1f",
        trans = "log10")
alk_fig("c2plus_fraction", "C2+\nfraction", "Figure_S7_c2plus_fraction", "%.2f", limits = c(0, 1))

# --- Figures S8-S11: bootstrap weight distributions, top 5 sources -------------
# Refit with keep_draws = TRUE (same seed and settings as script 01, so the
# draws are those summarized in the tables).
bt <- smma_bootstrap(res$input, B = cfg$B, seed = cfg$seed, kow = kow, sn = sn,
                     kow_args = cfg$kow_args, min_analytes = cfg$min_shared, keep_draws = TRUE)
stopifnot(isTRUE(all.equal(bt$per_source, res$boot$per_source)))
src_family <- function(s) family_of(catmap$category_maintext[match(s, catmap$source)])
labels <- c(`Factor 1` = "S8", `Factor 2` = "S9", `Factor 3` = "S10", `Factor 4` = "S11")
for (fid in names(labels)) {
  ms <- res$mix$sources[res$mix$sources$factor == fid, ]
  top <- utils::head(ms$source[order(-ms$weight)], 5)
  D <- bt$draws[[fid]][, top, drop = FALSE]; w <- bt$draw_weights[[fid]]
  d <- data.frame(source = rep(top, each = nrow(D)), value = as.vector(D), w = rep(w, length(top)))
  q <- do.call(rbind, lapply(top, function(s) {
    v <- smma:::weighted_quantile(D[, s], w)
    data.frame(source = s, median = v[1], lo = v[2], hi = v[3])
  }))
  lev <- rev(top)
  d$source <- factor(d$source, lev); q$source <- factor(q$source, lev)
  d$family <- factor(src_family(as.character(d$source)), family_levels)
  p <- ggplot(d, aes(x = .data$value, y = .data$source)) +
    geom_violin(aes(weight = .data$w, fill = .data$family), scale = "width", width = 0.8,
                colour = NA, alpha = 0.85, bw = 0.02) +
    geom_errorbar(data = q, aes(xmin = .data$lo, xmax = .data$hi, y = .data$source),
                  inherit.aes = FALSE, width = 0.2, linewidth = 0.4, orientation = "y") +
    geom_point(data = q, aes(x = .data$median, y = .data$source), inherit.aes = FALSE,
               size = 1.6, shape = 21, fill = "white", stroke = 0.6) +
    scale_x_continuous(limits = c(0, 1), labels = function(x) paste0(100 * x, "%")) +
    fill_family +
    labs(x = sprintf("SMMA weight across %s bootstrap draws (fit-weighted)",
                     format(cfg$B, big.mark = ",")),
         y = NULL, title = fid) +
    theme_pub + theme(legend.position = "bottom", plot.title = element_text(size = 9)) +
    guides(fill = guide_legend(nrow = 2))
  save_fig(p, sprintf("Figure_%s_bootstrap_weights_%s", labels[[fid]], sub(" ", "", fid)), 8)
}

# --- Figures S12-S15: observed versus reconstructed ----------------------------
labels <- c(`Factor 1` = "S12", `Factor 2` = "S13", `Factor 3` = "S14", `Factor 4` = "S15")
for (fid in names(labels)) {
  p <- plot_reconstruction(res$input, res$mix$sources, fid) +
    scale_fill_manual(values = c(Observed = "#52514e", Reconstructed = "#9ec5f4"), name = NULL) +
    theme(text = element_text(size = 8), legend.key.size = unit(3, "mm"),
          plot.title = element_text(size = 9), axis.text.x = element_text(size = 5, angle = 90,
                                                                      vjust = 0.5, hjust = 1))
  save_fig(p, sprintf("Figure_%s_reconstruction_%s", labels[[fid]], sub(" ", "", fid)), 8)
}

# --- Figures S24-S26: attributed sum of EPA16 by site class ---------------------
site_lev <- c("Urban", "Mixed", "Remote")
by_site <- function(tag) {
  t <- read.csv(file.path(cfg$out_dir, sprintf("attribution_mean_ci_by_site_class_%s.csv", tag)))
  t <- t[t$grouping == "EPA16" & t$mean > 0, ]
  t$group <- factor(t$group, site_lev)
  t
}
t17 <- by_site("17cat")
t17$family <- factor(family_of_detail[t17$category], family_levels)
t17$category <- factor(t17$category, rev(detail_levels))
pS24 <- ggplot(t17, aes(x = .data$share_pct, y = .data$category, fill = .data$family)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = pct_lab(.data$share_pct)), hjust = -0.15, size = 2.2, colour = "grey15") +
  facet_wrap(~ group, nrow = 1) +
  scale_x_continuous(limits = c(0, 45), expand = c(0, 0), labels = function(x) paste0(x, "%")) +
  scale_y_discrete(labels = dash) + fill_family +
  labs(x = expression("Share of mean attributed " * Sigma * "16 EPA PACs"), y = NULL) +
  theme_pub + theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 3))
save_fig(pS24, "Figure_S24_composition_17cat", 10)

conc_fig <- function(t, ylev) {
  t$category <- factor(t$category, rev(ylev))
  ggplot(t, aes(x = .data$mean, y = .data$category, fill = .data$family)) +
    geom_errorbar(aes(xmin = pmax(.data$lo, 0.001), xmax = .data$hi), width = 0.3,
                  linewidth = 0.35, colour = "grey35", orientation = "y") +
    geom_point(shape = 21, size = 2, stroke = 0.3, colour = "grey15") +
    facet_wrap(~ group, nrow = 1) +
    scale_x_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE, drop0trailing = TRUE)) +
    annotation_logticks(sides = "b", linewidth = 0.2, short = unit(0.5, "mm"),
                        mid = unit(1, "mm"), long = unit(1.5, "mm")) +
    scale_y_discrete(labels = dash) + fill_family +
    labs(x = expression("Mean attributed " * Sigma * "16 EPA PACs (ng g"^-1 * " dw) with 95% CI"),
         y = NULL) +
    theme_pub + theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 3))
}
save_fig(conc_fig(t17, detail_levels), "Figure_S25_concentration_17cat", 10)
t9 <- by_site("9cat")
t9$family <- factor(family_of(t9$category), family_levels)
save_fig(conc_fig(t9, family_levels[-8]), "Figure_S26_concentration_9cat", 7.5)

# --- Figure S21: PMF-then-SMMA versus SMMA applied directly to samples ---------
sites <- read.csv(file.path(cfg$out_dir, "site_classification.csv"))
pmf <- read.csv(file.path(cfg$out_dir, "attribution_median_share_by_site_class_9cat.csv"))
pmf <- pmf[pmf$grouping == "EPA16", ]
pt <- read.csv(file.path(cfg$out_dir, "smma_direct_samples_table_s10.csv"))
pc <- read.csv(file.path(cfg$out_dir, "smma_direct_samples_per_category.csv"))   # per sample
pc$site_class <- sites$site_class[match(pc$factor, sites$sample)]
bs <- stats::aggregate(cbind(median, p025, p975) ~ site_class + category, data = pc,
                       FUN = stats::median)
meth <- c("PMF then SMMA", "Direct SMMA, point estimate",
          "Direct SMMA, bootstrap median (95% CI)")
d21 <- rbind(
  data.frame(site_class = pmf$group, category = pmf$category, method = meth[1],
             value = pmf$median_pct, lo = NA, hi = NA),
  data.frame(site_class = pt$site_class, category = pt$category, method = meth[2],
             value = pt$median_pct, lo = NA, hi = NA),
  data.frame(site_class = bs$site_class, category = bs$category, method = meth[3],
             value = 100 * bs$median, lo = 100 * bs$p025, hi = 100 * bs$p975))
d21 <- d21[!grepl("^Light petroleum", d21$category), ]
d21$category <- factor(d21$category, rev(family_levels))
d21$site_class <- factor(d21$site_class, site_lev)
d21$method <- factor(d21$method, rev(meth))
pd <- position_dodge(width = 0.8)
pS21 <- ggplot(d21, aes(x = .data$value, y = .data$category, fill = .data$method)) +
  geom_col(position = pd, width = 0.75) +
  geom_errorbar(aes(xmin = .data$lo, xmax = .data$hi), position = pd, width = 0.3,
                linewidth = 0.3, colour = "grey20", orientation = "y", na.rm = TRUE) +
  facet_wrap(~ site_class, nrow = 1) +
  scale_fill_manual(values = stats::setNames(c("#0b0b0b", "#8a8984", "#d4d3cd"), meth),
                    breaks = meth, name = NULL) +
  scale_x_continuous(labels = function(x) paste0(x, "%")) +
  scale_y_discrete(labels = dash) +
  labs(x = "Median share of EPA16 burden", y = NULL) +
  theme_pub + theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 1))
save_fig(pS21, "Figure_S21_direct_vs_PMF_SMMA", 10)
# --- Figure S17: category share versus percent fines (Table S6) -------------------
s6 <- read.csv(file.path(cfg$out_dir, "share_vs_fines_table_s6.csv"))
s6$sig <- ifelse(s6$p_value >= 0.05, "Not significant",
                 ifelse(s6$r > 0, "Significant positive", "Significant negative"))
s6$category <- factor(s6$category, rev(c(family_levels[-8], "All petroleum")))
pS17 <- ggplot(s6, aes(x = .data$r, y = .data$category, fill = .data$sig)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_text(aes(label = sprintf("r = %.2f, p = %s", .data$r, signif(.data$p_value, 2)),
                hjust = ifelse(.data$r > 0, -0.08, 1.08)), size = 2.3, colour = "grey15") +
  scale_fill_manual(values = c(`Significant positive` = "#0ca30c",
                               `Significant negative` = "#d03b3b",
                               `Not significant` = "#b9b8b1"), name = NULL) +
  scale_x_continuous(limits = c(-0.6, 0.6)) + scale_y_discrete(labels = dash) +
  labs(x = sprintf("Pearson r, category share (EPA16) vs %% fines (n = %d)", s6$n[1]), y = NULL) +
  theme_pub + theme(legend.position = "bottom")
save_fig(pS17, "Figure_S17_share_vs_fines", 7.5)
cat("Figures written to", fig_dir, "\n")
