#' Plot SMMA category weights per factor
#'
#' @param boot An `smma_bootstrap` object, or the `categories` element of
#'   [smma_best_mixture()].
#' @param min_weight Hide categories whose median weight is below this value.
#' @return A ggplot object (requires \pkg{ggplot2}).
#' @export
plot_smma <- function(boot, min_weight = 0.005) {
  need_pkg("ggplot2", "for plotting")
  if (inherits(boot, "smma_bootstrap")) {
    d <- boot$per_category
    d$w <- d$median
  } else {
    d <- boot
    d$w <- d$weight
    d$p025 <- d$p975 <- NA_real_
  }
  d <- d[d$w >= min_weight, ]
  d$category <- stats::reorder(d$category, d$w)
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$w, y = .data$category)) +
    ggplot2::geom_col(fill = "#2a78d6", width = 0.7) +
    ggplot2::facet_wrap(~ factor, scales = "free_y") +
    ggplot2::scale_x_continuous(labels = function(x) paste0(round(100 * x), "%")) +
    ggplot2::labs(x = "SMMA weight", y = NULL) +
    ggplot2::theme_bw(base_size = 11)
  if (any(is.finite(d$p975))) {
    p <- p + ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$p025, xmax = .data$p975),
                                    width = 0.25, colour = "grey30", orientation = "y")
  }
  p
}

#' Plot observed versus reconstructed factor composition
#'
#' @param input An `smma_input` object.
#' @param mixture Data frame with `factor`, `source`, `weight`
#'   (e.g. `smma_best_mixture()$sources`).
#' @param factor_id Factor to plot.
#' @return A ggplot object.
#' @export
plot_reconstruction <- function(input, mixture, factor_id) {
  need_pkg("ggplot2", "for plotting")
  m <- mixture[mixture$factor == factor_id, ]
  w <- stats::setNames(numeric(ncol(input$L)), colnames(input$L))
  w[m$source] <- m$weight
  rec <- as.numeric(input$L %*% w)
  obs <- input$F[factor_id, ]
  rec <- rec / sum(rec[is.finite(obs)]); obs <- obs / sum(obs, na.rm = TRUE)
  d <- data.frame(analyte = rep(input$analytes, 2),
                  type = rep(c("Observed", "Reconstructed"), each = length(obs)),
                  fraction = c(obs, rec))
  d$analyte <- factor(d$analyte, levels = input$analytes)
  ggplot2::ggplot(d, ggplot2::aes(x = .data$analyte, y = .data$fraction, fill = .data$type)) +
    ggplot2::geom_col(position = "dodge", width = 0.8) +
    ggplot2::scale_fill_manual(values = c(Observed = "#52514e", Reconstructed = "#eb6834")) +
    ggplot2::labs(x = NULL, y = "Normalized fraction", fill = NULL, title = factor_id) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
                   legend.position = "top")
}
