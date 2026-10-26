# Heatmap of the criterion, or p-values for functional Box-Jenkins
#' @rdname select_order
#' @export
plot.farfima_order = function(x, main = NULL, ...) {
  tab = x$table
  method = attr(x, "method")

  if (method == "bj") {
    if (is.null(main)) {
      main = "Functional Box-Jenkins: p-values of the candidates tried"
    }
    graphics::barplot(tab$p_value, names.arg = paste0("(", tab$p, ",", tab$q, ")"), ylim = c(0, 1),
                      col = ifelse(tab$passed, "#1b9e77", "grey70"), border = NA, xlab = "(p, q)",
                      ylab = "p-value", main = main, ...)
    graphics::abline(h = attr(x, "alpha"), lty = 2)
    return(invisible(x))
  }

  score = if (method == "cv") tab[[if (attr(x, "cv_score") == "rmse") "RMSPE" else "MAPE"]] else tab$IC
  ps = seq(min(tab$p), max(tab$p))
  qs = seq(min(tab$q), max(tab$q))
  value = matrix(NA_real_, length(ps), length(qs))
  value[cbind(tab$p - min(ps) + 1, tab$q - min(qs) + 1)] = score
  ok = is.finite(value)
  colour_rank = value
  colour_rank[ok] = rank(value[ok])
  colour_rank[!ok] = NA
  if (is.null(main)) {
    main = switch(method, ic_full = "Functional IC", ic_auto = "Auto functional IC (candidates tried)",
                  cv = paste("CV", if (attr(x, "cv_score") == "rmse") "RMSPE" else "MAPE"))
  }

  graphics::image(ps, qs, colour_rank, col = grDevices::hcl.colors(max(sum(ok), 2), "YlOrRd", rev = TRUE),
                  xlab = "p", ylab = "q", main = main, axes = FALSE, ...)
  graphics::axis(1, at = ps)
  graphics::axis(2, at = qs)
  graphics::box()
  label = ifelse(ok, signif(value, 3), ifelse(is.na(value), "", "failed"))
  graphics::text(rep(ps, length(qs)), rep(qs, each = length(ps)), label, cex = 0.7)
  if (method == "ic_auto") {
    graphics::lines(x$search_path$p, x$search_path$q, type = "b", pch = 19, col = "grey20")
  }
  graphics::rect(x$p - 0.5, x$q - 0.5, x$p + 0.5, x$q + 0.5, lwd = 2.5)
  invisible(x)
}
