#' Attribute sample concentrations to library sources
#'
#' Propagates SMMA source weights back to samples. For sample *i*, source
#' *s* and analyte grouping *G*:
#'
#' \deqn{C(i,s,G) = \sum_f g(i,f) \, w(f,s) \sum_{a \in G} F(f,a)}
#'
#' where `g` are the receptor-model contributions (scores), `F` the factor
#' profiles in concentration units and `w` the SMMA weights (renormalized
#' within each factor over the retained sources).
#'
#' @param contributions Data frame of factor contributions, one row per
#'   sample (see [read_factor_contributions()]).
#' @param profiles Factor profiles in concentration units (see
#'   [read_factor_profiles()]); not the normalized `smma_input$F`.
#' @param mixture Source weights: the `sources` element of
#'   [smma_best_mixture()] or any data frame with `factor`, `source`,
#'   `weight` and a category column.
#' @param groupings Named list of analyte vectors, e.g. from
#'   [analyte_groupings()].
#' @param sample_col,factor_col Identifier columns in `contributions` and
#'   `profiles` (factor columns in `contributions` must match the factor
#'   identifiers in `profiles`).
#' @param category_col Category column in `mixture` used for aggregation;
#'   set `level = "source"` to keep individual sources.
#' @param category_map Optional data frame with columns `source` and
#'   `category` that replaces the categories of `mixture` (useful to
#'   regroup sources without refitting).
#' @param level `"category"` (default) or `"source"`.
#' @param negative How to treat negative contributions: `"drop"` (as in
#'   the accompanying analysis) or `"zero"`.
#' @return Long data frame with columns `sample`, `category` (or
#'   `source`), `grouping` and `value`. Every sample x category x grouping
#'   combination is present (zero when nothing was attributed), so means
#'   and medians are not biased by missing rows.
#' @examples
#' ex <- smma_example()
#' bt <- smma_bootstrap(smma_prepare(ex$profiles, ex$library), B = 30, seed = 1)
#' mix <- smma_best_mixture(bt)
#' att <- attribute_sources(ex$contributions, ex$profiles, mix$sources,
#'                          analyte_groupings(names(ex$profiles)))
#' head(att)
#' @export
attribute_sources <- function(contributions, profiles, mixture, groupings,
                              sample_col = "sample", factor_col = "factor",
                              category_col = "category", category_map = NULL,
                              level = c("category", "source"),
                              negative = c("drop", "zero")) {
  level <- match.arg(level)
  negative <- match.arg(negative)
  check_cols(mixture, c("factor", "source", "weight"), "mixture")
  check_cols(contributions, sample_col, "contributions")
  check_cols(profiles, factor_col, "profiles")
  if (!is.list(groupings) || is.null(names(groupings))) stop("groupings must be a named list.")

  mix <- mixture[is.finite(mixture$weight) & mixture$weight > 0 & !is.na(mixture$source), ]
  if (!is.null(category_map)) {
    check_cols(category_map, c("source", "category"), "category_map")
    mix$.cat <- category_map$category[match(trimws(mix$source), trimws(category_map$source))]
    if (anyNA(mix$.cat)) {
      warning("Sources without a category in category_map: ",
              paste(unique(mix$source[is.na(mix$.cat)]), collapse = ", "), call. = FALSE)
    }
  } else {
    check_cols(mixture, category_col, "mixture")
    mix$.cat <- mix[[category_col]]
  }
  mix$weight <- stats::ave(mix$weight, mix$factor, FUN = function(w) w / sum(w))
  key <- if (level == "category") mix$.cat else mix$source

  factors <- intersect(unique(mix$factor), setdiff(names(contributions), sample_col))
  if (!length(factors)) stop("No factor columns of 'contributions' match the factors in 'mixture'.")
  prof <- profiles[match(factors, as.character(profiles[[factor_col]])), , drop = FALSE]
  if (anyNA(prof[[factor_col]])) stop("Some factors are missing from 'profiles'.")

  # factor x grouping concentration sums
  gsum <- sapply(groupings, function(a) {
    a <- intersect(a, names(prof))
    if (!length(a)) return(rep(NA_real_, nrow(prof)))
    m <- as.matrix(as_numeric_cols(prof, a)[, a, drop = FALSE])
    rowSums(m, na.rm = TRUE)
  })
  gsum <- matrix(gsum, nrow = length(factors), dimnames = list(factors, names(groupings)))

  G <- as.matrix(as_numeric_cols(contributions, factors)[, factors, drop = FALSE])
  if (negative == "zero") G[G < 0] <- 0 else G[G < 0] <- NA_real_
  samples <- as.character(contributions[[sample_col]])
  levels_key <- sort(unique(key[!is.na(key)]), method = "radix")

  # weight matrix factors x keys
  Wk <- matrix(0, length(factors), length(levels_key), dimnames = list(factors, levels_key))
  for (j in seq_len(nrow(mix))) {
    if (mix$factor[j] %in% factors && !is.na(key[j])) {
      Wk[mix$factor[j], key[j]] <- Wk[mix$factor[j], key[j]] + mix$weight[j]
    }
  }
  out <- lapply(names(groupings), function(gname) {
    # sample x key = sum_f G[i,f] * gsum[f] * Wk[f,k]; NA contributions skipped
    Gf <- G
    Gf[is.na(Gf)] <- 0
    val <- Gf %*% (gsum[, gname] * Wk)
    data.frame(sample = rep(samples, times = ncol(val)),
               key = rep(colnames(val), each = nrow(val)),
               grouping = gname, value = as.vector(val), stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, out)
  names(res)[2] <- level
  res
}

#' Summarise attributed concentrations by group
#'
#' Computes, for each group (e.g. urban/mixed/remote site class) x
#' category x analyte grouping, the number of samples, the mean
#' attributed concentration with a percentile-bootstrap confidence
#' interval (resampling samples), the median, and the share of the
#' group's mean burden.
#'
#' @param attribution Output of [attribute_sources()].
#' @param groups Data frame mapping `sample` to a grouping column, or
#'   `NULL` to summarise all samples together.
#' @param group_col Name of the grouping column in `groups`.
#' @param level Column of `attribution` holding categories (or sources).
#' @param B Bootstrap resamples for the confidence interval (0 to skip).
#' @param conf Confidence level.
#' @param seed Optional seed.
#' @return Data frame with `group`, `level`, `grouping`, `n`, `mean`,
#'   `lo`, `hi`, `median`, `share_pct`.
#' @export
summarise_attribution <- function(attribution, groups = NULL, group_col = "group",
                                  level = "category", B = 2000, conf = 0.95, seed = NULL) {
  att <- attribution
  if (is.null(groups)) {
    att$.group <- "All"
  } else {
    check_cols(groups, c("sample", group_col), "groups")
    att$.group <- groups[[group_col]][match(att$sample, groups$sample)]
    att <- att[!is.na(att$.group), ]
  }
  if (!is.null(seed)) set.seed(seed)
  a <- (1 - conf) / 2
  keys <- unique(att[, c(".group", level, "grouping")])
  keys <- keys[order(as.character(keys$.group), keys[[level]], keys$grouping, method = "radix"), ]
  rows <- lapply(seq_len(nrow(keys)), function(k) {
    x <- att$value[att$.group == keys$.group[k] & att[[level]] == keys[[level]][k] &
                     att$grouping == keys$grouping[k]]
    ci <- if (B > 0 && length(x) > 1) {
      bm <- replicate(B, mean(sample(x, length(x), replace = TRUE)))
      stats::quantile(bm, c(a, 1 - a), names = FALSE)
    } else c(NA_real_, NA_real_)
    data.frame(group = keys$.group[k], key = keys[[level]][k], grouping = keys$grouping[k],
               n = length(x), mean = mean(x), lo = ci[1], hi = ci[2],
               median = stats::median(x), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  tot <- stats::ave(out$mean, out$group, out$grouping, FUN = function(v) sum(v, na.rm = TRUE))
  out$share_pct <- 100 * out$mean / tot
  names(out)[names(out) == "key"] <- level
  rownames(out) <- NULL
  out
}

#' Per-sample percentage shares and their medians
#'
#' Expresses each sample's attributed burden as percentages by category
#' (within each analyte grouping) and returns the median percentage per
#' group x category x grouping.
#'
#' @inheritParams summarise_attribution
#' @return Data frame with `group`, `level`, `grouping`, `n`, `median_pct`.
#' @export
attribution_shares <- function(attribution, groups = NULL, group_col = "group",
                               level = "category") {
  att <- attribution
  tot <- stats::ave(att$value, att$sample, att$grouping, FUN = function(v) sum(v, na.rm = TRUE))
  att$pct <- ifelse(tot > 0, 100 * att$value / tot, NA_real_)
  att$.group <- if (is.null(groups)) "All" else groups[[group_col]][match(att$sample, groups$sample)]
  att <- att[!is.na(att$.group), ]
  agg <- stats::aggregate(att$pct, by = list(group = att$.group, key = att[[level]],
                                             grouping = att$grouping),
                          FUN = function(v) c(n = sum(!is.na(v)), median = stats::median(v, na.rm = TRUE)))
  out <- data.frame(agg[, 1:3], n = agg$x[, "n"], median_pct = agg$x[, "median"])
  names(out)[2] <- level
  out[order(out$group, out$grouping, -out$median_pct), ]
}
