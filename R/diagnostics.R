#' Classic PAC diagnostic ratios
#'
#' Computes Ant/(Ant+Phe), Fla/(Fla+Pyr), BaA/(BaA+Chr) and
#' IcdP/(IcdP+BghiP) for each row of a profile or sample table, using the
#' compound columns defined in the analyte scheme.
#'
#' @param x Data frame of factor profiles or samples.
#' @param id_col Identifier column.
#' @param scheme Analyte scheme ([default_analyte_scheme()]).
#' @return Data frame with `id_col` and one column per ratio.
#' @examples
#' prof <- read_factor_profiles(system.file("extdata", "bc_factor_profiles.csv",
#'                                          package = "smma"))
#' diagnostic_ratios(prof)
#' @export
diagnostic_ratios <- function(x, id_col = "factor", scheme = default_analyte_scheme()) {
  r <- scheme$ratio
  need <- c("ant", "phe", "fla", "pyr", "baa", "chr", "icdp", "bghip")
  if (!all(need %in% names(r))) stop("scheme$ratio must define: ", paste(need, collapse = ", "))
  check_cols(x, c(id_col, unname(r[need])), "x")
  v <- function(k) suppressWarnings(as.numeric(x[[r[[k]]]]))
  frac <- function(a, b) ifelse(is.finite(a + b) & (a + b) > 0, a / (a + b), NA_real_)
  out <- data.frame(x[[id_col]],
                    frac(v("ant"), v("phe")), frac(v("fla"), v("pyr")),
                    frac(v("baa"), v("chr")), frac(v("icdp"), v("bghip")),
                    stringsAsFactors = FALSE)
  names(out) <- c(id_col, "Ant_AntPhe", "Fla_FlaPyr", "BaA_BaAChr", "IcdP_IcdPBghiP")
  out
}

#' Alkylation indices within homologue series
#'
#' For every homologue series in the analyte scheme, returns the parent
#' fraction (C0 / sum), the alkylated-to-parent ratio ((C1+...)/C0) and
#' the C2+ fraction ((C2+...)/sum). High parent fractions indicate
#' pyrogenic input; high alkyl/parent ratios and C2+ fractions indicate
#' petrogenic input.
#'
#' @inheritParams diagnostic_ratios
#' @return Long data frame with `id_col`, `series`, `parent_fraction`,
#'   `alkyl_parent_ratio`, `c2plus_fraction`.
#' @export
alkylation_indices <- function(x, id_col = "factor", scheme = default_analyte_scheme()) {
  check_cols(x, id_col, "x")
  rows <- lapply(names(scheme$homologues), function(s) {
    cols <- intersect(scheme$homologues[[s]], names(x))
    if (length(cols) < 2) return(NULL)
    m <- as.matrix(as_numeric_cols(x, cols)[, cols, drop = FALSE])
    tot <- rowSums(m, na.rm = TRUE); p <- m[, 1]
    alk <- rowSums(m[, -1, drop = FALSE], na.rm = TRUE)
    c2 <- if (ncol(m) >= 3) rowSums(m[, -(1:2), drop = FALSE], na.rm = TRUE) else NA_real_
    data.frame(x[[id_col]], s,
               ifelse(tot > 0, p / tot, NA_real_),
               ifelse(p > 0, alk / p, NA_real_),
               ifelse(tot > 0, c2 / tot, NA_real_), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  names(out) <- c(id_col, "series", "parent_fraction", "alkyl_parent_ratio", "c2plus_fraction")
  out
}

#' Distances between factor profiles and library profiles
#'
#' Compares every factor with every library source over the analytes
#' they share (both renormalized to unit sum on that subset) using cosine
#' similarity, correlation distance (1 - Pearson r), Bray-Curtis
#' distance, Euclidean distance and mean chi-square. Optionally adds a
#' weighted chi-square on log-transformed, total-scaled profiles with
#' per-analyte weights (e.g. Kow x S/N).
#'
#' @param profiles Factor profiles (one row per factor).
#' @param source_library Source library.
#' @param factor_col,source_col,category_col Identifier columns.
#' @param analytes Analytes to use (default: shared columns matching `analyte_pattern`).
#' @param analyte_pattern Pattern used when `analytes` is `NULL`.
#' @param min_shared Minimum shared analytes for a comparison.
#' @param weights Optional function `function(factor_id, analytes)`
#'   returning per-analyte weights for the weighted chi-square.
#' @return Data frame, one row per factor x source, with `n_shared` and
#'   the metrics. Use [best_matches()] to rank.
#' @export
profile_distances <- function(profiles, source_library, factor_col = "factor",
                              source_col = "source", category_col = "category",
                              analytes = NULL, analyte_pattern = "^(pah_|tph_|biom_)",
                              min_shared = 12, weights = NULL) {
  if (is.null(analytes)) {
    analytes <- intersect(grep(analyte_pattern, names(profiles), value = TRUE),
                          grep(analyte_pattern, names(source_library), value = TRUE))
  }
  P <- as.matrix(as_numeric_cols(profiles, analytes)[, analytes, drop = FALSE])
  S <- as.matrix(as_numeric_cols(source_library, analytes)[, analytes, drop = FALSE])
  res <- list(); k <- 0
  for (i in seq_len(nrow(P))) {
    fid <- as.character(profiles[[factor_col]][i])
    wv <- if (is.null(weights)) NULL else weights(fid, analytes)
    for (j in seq_len(nrow(S))) {
      ok <- is.finite(P[i, ]) & is.finite(S[j, ])
      n <- sum(ok)
      met <- c(cosine = NA, corr_dist = NA, bray_curtis = NA, euclidean = NA,
               chisq_mean = NA, chisq_weighted = NA)
      if (n >= min_shared && sum(P[i, ok]) > 0 && sum(S[j, ok]) > 0) {
        f <- P[i, ok] / sum(P[i, ok]); s <- S[j, ok] / sum(S[j, ok])
        den <- sqrt(sum(f^2)) * sqrt(sum(s^2))
        met["cosine"] <- if (den > 0) sum(f * s) / den else NA
        cr <- suppressWarnings(stats::cor(f, s))
        met["corr_dist"] <- if (is.finite(cr)) 1 - cr else NA
        met["bray_curtis"] <- sum(abs(f - s)) / sum(f + s)
        met["euclidean"] <- sqrt(sum((f - s)^2))
        met["chisq_mean"] <- if (any(s > 0)) mean((f[s > 0] - s[s > 0])^2 / s[s > 0]) else NA
        if (!is.null(wv)) {
          w <- wv[ok]; w[!is.finite(w) | w <= 0] <- 1
          raw <- P[i, ok]; e <- s; e[e == 0] <- 1e-4
          o <- log(raw + 1); ex <- log(e * sum(raw) + 1)
          met["chisq_weighted"] <- sum(w * (o - ex)^2 / pmax(ex, .Machine$double.eps)) / length(o)
        }
      }
      k <- k + 1
      res[[k]] <- data.frame(factor = fid, source = as.character(source_library[[source_col]][j]),
                             category = as.character(source_library[[category_col]][j]),
                             n_shared = n, t(met), stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, res)
}

#' Best-matching library sources per factor and metric
#'
#' @param distances Output of [profile_distances()].
#' @param n Number of matches to return per factor and metric.
#' @return Long data frame: `factor`, `metric`, `rank`, `source`,
#'   `category`, `score`, `n_shared`.
#' @export
best_matches <- function(distances, n = 1) {
  metrics <- intersect(c("cosine", "corr_dist", "bray_curtis", "euclidean", "chisq_mean",
                         "chisq_weighted"), names(distances))
  out <- lapply(metrics, function(m) {
    d <- distances[is.finite(distances[[m]]), ]
    do.call(rbind, lapply(split(d, d$factor), function(df) {
      o <- if (m == "cosine") order(-df[[m]], -df$n_shared) else order(df[[m]], -df$n_shared)
      df <- df[o[seq_len(min(n, nrow(df)))], ]
      data.frame(factor = df$factor, metric = m, rank = seq_len(nrow(df)), source = df$source,
                 category = df$category, score = df[[m]], n_shared = df$n_shared,
                 stringsAsFactors = FALSE)
    }))
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}
