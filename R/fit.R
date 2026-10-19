# Core weighted NNLS on unit-sum profiles. Returns mixing coefficients
# renormalized to sum to 1.
nnls_unit <- function(f, L, w = NULL) {
  if (!is.null(w)) {
    ok <- is.finite(w) & w > 0
    if (!any(ok)) w <- rep(1, length(f))
    w[!ok] <- 0
    if (sum(w) > 0) {
      w <- w / mean(w[w > 0])
      sw <- sqrt(w)
      f <- f * sw
      L <- L * sw
    }
  }
  s <- stats::coef(nnls::nnls(L, f))
  s <- as.numeric(s)
  if (sum(s) > 0) s <- s / sum(s)
  s
}

# Restrict to analytes that are finite in the factor and non-zero in at
# least one source, then renormalize factor and sources on that subset.
fit_subset <- function(f, L, idx, w) {
  fv <- f[idx]; Lv <- L[idx, , drop = FALSE]
  fv[!is.finite(fv)] <- 0
  Lv[!is.finite(Lv)] <- 0
  sf <- sum(fv)
  if (sf > 0) fv <- fv / sf
  cs <- colSums(Lv)
  cs[cs <= 0] <- NA_real_
  Lv <- sweep(Lv, 2, cs, "/")
  Lv[!is.finite(Lv)] <- 0
  s <- nnls_unit(fv, Lv, w)
  recon <- as.numeric(Lv %*% s)
  list(s = s, rel_rmse = rel_rmse(fv - recon, fv, w), f = fv, recon = recon)
}

factor_weights <- function(input, factor_id, kow = NULL, sn = NULL,
                           kow_args = list(), sn_factor_col = "factor") {
  a <- rownames(input$L)
  kw <- do.call(kow_weights, c(list(analytes = a, kow = kow), kow_args))
  sw <- sn_weights(sn, factor_id, a, factor_col = sn_factor_col)
  combine_weights(kw, sw)
}

valid_analytes <- function(input, factor_id) {
  f <- input$F[factor_id, ]
  which(is.finite(f) & rowSums(input$L) > 0)
}

#' Fit each factor as a mixture of library sources (single fit)
#'
#' Solves, for every factor, the weighted non-negative least-squares
#' problem `min || sqrt(W) (f - L s) ||^2, s >= 0` on the analytes shared by
#' the factor and the library, then renormalizes `s` to sum to 1. This is
#' the point estimate; use [smma_bootstrap()] for uncertainty.
#'
#' @param input An `smma_input` object from [smma_prepare()].
#' @param kow Optional named Kow vector (see [read_kow()]).
#' @param sn Optional factor-specific S/N table (see [read_sn()]).
#' @param kow_args List of further arguments to [kow_weights()]
#'   (e.g. `list(input_scale = "log10", matrix = "sediment")`).
#' @param sn_factor_col Factor column name in `sn`.
#' @param min_analytes Minimum number of usable analytes per factor.
#' @return Data frame with columns `factor`, `source`, `category`,
#'   `weight` and `rel_rmse` (relative RMSE of the reconstruction).
#' @examples
#' ex <- smma_example()
#' fit <- smma_fit(smma_prepare(ex$profiles, ex$library))
#' head(fit[order(-fit$weight), ])
#' @export
smma_fit <- function(input, kow = NULL, sn = NULL, kow_args = list(),
                     sn_factor_col = "factor", min_analytes = 12) {
  stopifnot(inherits(input, "smma_input"))
  out <- lapply(rownames(input$F), function(fid) {
    idx <- valid_analytes(input, fid)
    if (length(idx) < min_analytes) {
      warning("Factor ", fid, ": too few usable analytes; skipped.", call. = FALSE)
      return(NULL)
    }
    w <- factor_weights(input, fid, kow, sn, kow_args, sn_factor_col)[idx]
    r <- fit_subset(input$F[fid, ], input$L, idx, w)
    data.frame(factor = fid, source = colnames(input$L), weight = r$s,
               rel_rmse = r$rel_rmse, stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, out)
  res$category <- input$categories$category[match(res$source, input$categories$source)]
  res[, c("factor", "source", "category", "weight", "rel_rmse")]
}

#' Forward selection of a small number of sources
#'
#' Greedy selection of up to `max_k` library sources per factor, adding at
#' each step the source that most reduces the relative RMSE, and stopping
#' once the relative RMSE falls below `threshold`.
#'
#' @inheritParams smma_fit
#' @param max_k Maximum number of sources.
#' @param threshold Relative RMSE at which selection stops.
#' @return Data frame with `factor`, `step`, `source`, `category`,
#'   `weight` (final weights of the selected set) and `rel_rmse` (after
#'   each step).
#' @export
smma_forward_select <- function(input, max_k = 3, threshold = 0.10, min_analytes = 12) {
  stopifnot(inherits(input, "smma_input"))
  out <- lapply(rownames(input$F), function(fid) {
    idx <- valid_analytes(input, fid)
    if (length(idx) < min_analytes) return(NULL)
    fv <- input$F[fid, idx]; Lv <- input$L[idx, , drop = FALSE]
    remaining <- colnames(Lv); sel <- character(); path <- numeric(); best_w <- NULL
    for (k in seq_len(max_k)) {
      cand <- lapply(remaining, function(sid) {
        cols <- c(sel, sid)
        w <- stats::coef(nnls::nnls(Lv[, cols, drop = FALSE], fv))
        list(sid = sid, w = w,
             rel = rel_rmse(fv - as.numeric(Lv[, cols, drop = FALSE] %*% w), fv))
      })
      b <- cand[[which.min(vapply(cand, `[[`, 0, "rel"))]]
      sel <- c(sel, b$sid); path <- c(path, b$rel); best_w <- b$w
      remaining <- setdiff(remaining, b$sid)
      if (b$rel <= threshold || !length(remaining)) break
    }
    if (sum(best_w) > 0) best_w <- best_w / sum(best_w)
    data.frame(factor = fid, step = seq_along(sel), source = sel, weight = best_w,
               rel_rmse = path, stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, out)
  res$category <- input$categories$category[match(res$source, input$categories$source)]
  res[, c("factor", "step", "source", "category", "weight", "rel_rmse")]
}
