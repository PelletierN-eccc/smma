#' Bootstrapped Source Mixture Matching Analysis
#'
#' For each factor, repeatedly resamples the shared analytes with
#' replacement (same number of analytes, uniform probability),
#' renormalizes the factor and every library profile on the resampled
#' set, and refits the weighted NNLS mixture. Draws are weighted by the
#' inverse of their weighted relative RMSE, and mixing coefficients are
#' summarized per source and per category as weighted medians, 95%
#' weighted percentile intervals, weighted means and selection rates
#' (weighted share of draws in which the source received a positive
#' coefficient).
#'
#' Per-analyte weights are the product of [kow_weights()] and the
#' factor-specific S/N from `sn` ([sn_weights()]); non-finite or
#' non-positive products are set to 1 (neutral).
#'
#' Random numbers: when `seed` is given, one seed per factor is derived
#' from it, so results for a factor do not depend on how many other
#' factors are fitted or on their order.
#'
#' @inheritParams smma_fit
#' @param B Number of bootstrap draws per factor.
#' @param seed Integer seed, or `NULL` to use the current RNG state.
#' @param keep_best Optional proportion in (0, 1]: keep only this share of
#'   draws with the lowest relative RMSE (at least `keep_min` draws).
#' @param keep_min Minimum number of draws kept when `keep_best` is set.
#' @param keep_draws Store the full matrix of draws and their weights
#'   (needed for [smma_group_share()] and distribution plots; can be large).
#' @param progress Show a progress bar.
#' @return An object of class `smma_bootstrap`, a list with data frames
#'   `per_source` and `per_category` (columns `factor`, `source` or
#'   `category`, `median`, `p025`, `p975`, `mean`, `sel_rate`), `fit`
#'   (relative RMSE summary per factor), optionally `draws`, and `settings`.
#' @examples
#' ex <- smma_example()
#' inp <- smma_prepare(ex$profiles, ex$library)
#' bt <- smma_bootstrap(inp, B = 50, seed = 1)
#' head(bt$per_category)
#' @export
smma_bootstrap <- function(input, B = 1000, seed = NULL, kow = NULL, sn = NULL,
                           kow_args = list(), sn_factor_col = "factor",
                           min_analytes = 12, keep_best = NULL, keep_min = 20,
                           keep_draws = FALSE, progress = interactive()) {
  stopifnot(inherits(input, "smma_input"))
  fac_ids <- rownames(input$F)
  seeds <- NULL
  if (!is.null(seed)) {
    set.seed(seed)
    seeds <- stats::setNames(sample.int(.Machine$integer.max, length(fac_ids)), fac_ids)
  }
  cats <- input$categories
  src2cat <- stats::setNames(cats$category, cats$source)

  res <- lapply(fac_ids, function(fid) {
    idx <- valid_analytes(input, fid)
    if (length(idx) < min_analytes) {
      warning("Factor ", fid, ": too few usable analytes; skipped.", call. = FALSE)
      return(NULL)
    }
    if (!is.null(seeds)) set.seed(seeds[[fid]])
    w_all <- factor_weights(input, fid, kow, sn, kow_args, sn_factor_col)
    f <- input$F[fid, ]
    D <- matrix(NA_real_, B, ncol(input$L), dimnames = list(NULL, colnames(input$L)))
    rr <- rep(NA_real_, B)
    pb <- if (progress) utils::txtProgressBar(0, B, style = 3) else NULL
    for (b in seq_len(B)) {
      ib <- sample(idx, length(idx), replace = TRUE)
      r <- fit_subset(f, input$L, ib, w_all[ib])
      D[b, ] <- r$s
      rr[b] <- r$rel_rmse
      if (!is.null(pb) && b %% 50 == 0) utils::setTxtProgressBar(pb, b)
    }
    if (!is.null(pb)) close(pb)

    wf <- 1 / (rr + 1e-6)
    wf[!is.finite(wf)] <- 0
    if (sum(wf) <= 0) wf[] <- 1
    wf <- wf / sum(wf)
    keep <- seq_len(B)
    if (!is.null(keep_best)) {
      if (!is.finite(keep_best) || keep_best <= 0 || keep_best > 1)
        stop("keep_best must be in (0, 1].", call. = FALSE)
      n <- min(B, max(keep_min, ceiling(B * keep_best), 1))
      keep <- order(rr)[seq_len(n)]
    }
    Du <- D[keep, , drop = FALSE]
    wu <- wf[keep]
    if (sum(wu) <= 0) wu[] <- 1
    wu <- wu / sum(wu)

    summarise_draws <- function(M, id_name) {
      rows <- lapply(colnames(M), function(nm) {
        x <- M[, nm]
        w <- wu; w[!is.finite(x)] <- 0
        if (sum(w) <= 0) w[is.finite(x)] <- 1
        q <- weighted_quantile(x, w)
        data.frame(nm, q[1], q[2], q[3], weighted_mean(x, w),
                   weighted_mean(as.numeric(x > 0), w), stringsAsFactors = FALSE)
      })
      out <- do.call(rbind, rows)
      names(out) <- c(id_name, "median", "p025", "p975", "mean", "sel_rate")
      out
    }
    ps <- summarise_draws(Du, "source")
    ps$category <- unname(src2cat[ps$source])
    Mc <- t(apply(Du, 1, function(w) tapply(w, src2cat[colnames(Du)], sum, na.rm = TRUE)))
    if (length(unique(src2cat)) == 1) {
      Mc <- matrix(Mc, ncol = 1, dimnames = list(NULL, unique(src2cat)))
    }
    pc <- summarise_draws(Mc[, sort(colnames(Mc)), drop = FALSE], "category")
    list(per_source = cbind(factor = fid, ps, stringsAsFactors = FALSE),
         per_category = cbind(factor = fid, pc, stringsAsFactors = FALSE),
         fit = data.frame(factor = fid, n_analytes = length(idx), n_draws = length(keep),
                          rel_rmse_median = stats::median(rr[keep], na.rm = TRUE),
                          rel_rmse_p05 = stats::quantile(rr[keep], 0.05, na.rm = TRUE, names = FALSE),
                          rel_rmse_p95 = stats::quantile(rr[keep], 0.95, na.rm = TRUE, names = FALSE)),
         draws = if (keep_draws) Du else NULL,
         draw_weights = if (keep_draws) wu else NULL)
  })
  names(res) <- fac_ids
  res <- res[!vapply(res, is.null, logical(1))]
  rbind_el <- function(el) {
    out <- do.call(rbind, unname(lapply(res, `[[`, el)))
    rownames(out) <- NULL
    out
  }
  structure(list(
    per_source = rbind_el("per_source"),
    per_category = rbind_el("per_category"),
    fit = rbind_el("fit"),
    draws = if (keep_draws) lapply(res, `[[`, "draws") else NULL,
    draw_weights = if (keep_draws) lapply(res, `[[`, "draw_weights") else NULL,
    settings = list(B = B, seed = seed, keep_best = keep_best, min_analytes = min_analytes)
  ), class = "smma_bootstrap")
}

#' Print an smma_bootstrap object
#'
#' @param x An `smma_bootstrap` object.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.smma_bootstrap <- function(x, ...) {
  cat(sprintf("<smma_bootstrap> B = %d per factor, %d factors\n", x$settings$B, nrow(x$fit)))
  top <- x$per_category[x$per_category$median > 0.05, c("factor", "category", "median", "p025", "p975")]
  top <- top[order(top$factor, -top$median), ]
  top[, 3:5] <- round(top[, 3:5], 3)
  print(top, row.names = FALSE)
  invisible(x)
}

#' Best-guess source mixture from bootstrap results
#'
#' Keeps, within each factor, the sources whose selection rate is at
#' least `min_sel_rate` or whose median (or mean) weight is at least
#' `min_weight` (optionally only the `top_n` largest), and renormalizes
#' their weights to sum to 1. The top-ranked source is always kept.
#'
#' @param boot An `smma_bootstrap` object.
#' @param min_sel_rate,min_weight Retention thresholds.
#' @param top_n Optional maximum number of sources per factor.
#' @param stat Statistic used as the weight: `"median"` or `"mean"`.
#' @return List with `sources` (factor, source, category, weight and the
#'   bootstrap summary columns) and `categories` (weight summed by category).
#' @examples
#' ex <- smma_example()
#' bt <- smma_bootstrap(smma_prepare(ex$profiles, ex$library), B = 50, seed = 1)
#' smma_best_mixture(bt)$categories
#' @export
smma_best_mixture <- function(boot, min_sel_rate = 0.20, min_weight = 0.01,
                              top_n = NULL, stat = c("median", "mean")) {
  stat <- match.arg(stat)
  ps <- if (inherits(boot, "smma_bootstrap")) boot$per_source else boot
  ps$base <- ps[[stat]]
  ps$base[!is.finite(ps$base)] <- 0
  ps$sel_rate[!is.finite(ps$sel_rate)] <- 0
  out <- do.call(rbind, lapply(split(ps, ps$factor), function(d) {
    d <- d[order(-d$base), ]
    d$rank <- seq_len(nrow(d))
    keep <- (d$sel_rate >= min_sel_rate | d$base >= min_weight) &
      (if (is.null(top_n)) TRUE else d$rank <= top_n)
    d <- d[keep | d$rank == 1, ]
    d$weight <- if (sum(d$base) > 0) d$base / sum(d$base) else 0
    d
  }))
  out <- out[, c("factor", "source", "category", "weight", "median", "p025", "p975", "mean", "sel_rate")]
  rownames(out) <- NULL
  cats <- stats::aggregate(weight ~ factor + category, data = out, FUN = sum)
  cats <- cats[order(cats$factor, -cats$weight), ]
  rownames(cats) <- NULL
  list(sources = out, categories = cats)
}

#' Bootstrap share of a group of sources
#'
#' Sums, within every bootstrap draw, the weights of the sources assigned
#' to each group (e.g. petrogenic versus other) and summarises the summed
#' weight with the same fit-weighted statistics as [smma_bootstrap()].
#'
#' @param boot An `smma_bootstrap` object created with `keep_draws = TRUE`.
#' @param groups Named character vector mapping source identifiers to groups,
#'   or a function of the category returning the group.
#' @param categories Data frame `source`, `category` (e.g. `input$categories`),
#'   required when `groups` is a function.
#' @return Data frame `factor`, `group`, `median`, `p025`, `p975`, `mean`.
#' @export
smma_group_share <- function(boot, groups, categories = NULL) {
  if (is.null(boot$draws)) stop("Run smma_bootstrap(..., keep_draws = TRUE) first.", call. = FALSE)
  out <- lapply(names(boot$draws), function(fid) {
    D <- boot$draws[[fid]]; w <- boot$draw_weights[[fid]]
    g <- if (is.function(groups)) {
      if (is.null(categories)) stop("'categories' is required when 'groups' is a function.")
      groups(categories$category[match(colnames(D), categories$source)])
    } else groups[colnames(D)]
    do.call(rbind, lapply(sort(unique(stats::na.omit(g))), function(gr) {
      x <- rowSums(D[, which(g == gr), drop = FALSE])
      q <- weighted_quantile(x, w)
      data.frame(factor = fid, group = gr, median = q[1], p025 = q[2], p975 = q[3],
                 mean = weighted_mean(x, w), stringsAsFactors = FALSE)
    }))
  })
  do.call(rbind, out)
}
