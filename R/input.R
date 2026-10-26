# Turn a matrix, data frame or fts object into a grid-by-time numeric matrix
as_curves = function(x, name, n_grid = NULL) {
  if (inherits(x, "fts")) {
    y = as.matrix(x$y)
    rownames(y) = x$x
    x = y
  }
  if (is.data.frame(x)) {
    x = as.matrix(x)
  }

  if (!is.matrix(x) || !is.numeric(x)) {
    stop(sprintf("`%s` must be a numeric matrix (grid points by time) or an `fts` object.", name), call. = FALSE)
  }
  if (!is.null(n_grid) && nrow(x) != n_grid) {
    stop(sprintf("`%s` has %d grid points (rows) but the training curves have %d.", name, nrow(x), n_grid),
         call. = FALSE)
  }
  bad = sum(!is.finite(x))
  if (bad > 0) {
    stop(sprintf("`%s` has %d missing or infinite values.", name, bad), call. = FALSE)
  }

  x
}

# Stop when trimming leaves too few training curves
check_trim = function(n_curves, trim, d) {
  if (d != 0 && n_curves - trim < 10) {
    stop(sprintf("trim = %d leaves %d of the %d training curves; at least 10 are needed. Are the curves grid points by time?",
                 trim, max(0, n_curves - trim), n_curves), call. = FALSE)
  }
}

# Warn when a sample does not start right after the previous one
check_follow = function(prev, nxt, prev_name, next_name) {
  a = suppressWarnings(as.numeric(colnames(prev)))
  b = suppressWarnings(as.numeric(colnames(nxt)))
  if (length(a) < 2 || length(b) < 1 || any(!is.finite(a)) || any(!is.finite(b))) {
    return(invisible(NULL))
  }
  step = a[length(a)] - a[length(a) - 1]
  if (abs(b[1] - (a[length(a)] + step)) > 1e-8) {
    warning(sprintf("`%s` starts at %s rather than %s. The samples must follow each other in time.",
                    next_name, colnames(nxt)[1], format(a[length(a)] + step)), call. = FALSE)
  }
  invisible(NULL)
}

# Stop on unknown entries in memory_args or order_args
check_args = function(args, fun, name, exclude) {
  allowed = setdiff(names(formals(fun)), exclude)
  unknown = setdiff(names(args), allowed)
  if (length(unknown) > 0 || (length(args) > 0 && is.null(names(args)))) {
    stop(sprintf("Unknown entries in `%s`: %s. Allowed: %s.", name,
                 paste(if (length(unknown)) unknown else "<unnamed>", collapse = ", "),
                 paste(allowed, collapse = ", ")), call. = FALSE)
  }
}

# Time labels of the curves, or their positions when the columns have no names
origin_labels = function(X) {
  if (!is.null(colnames(X)) && all(nzchar(colnames(X)))) colnames(X) else as.character(seq_len(ncol(X)))
}

# Text progress bar, shown only in interactive sessions
progress_bar = function(verbose, n) {
  if (verbose && interactive() && n > 1) utils::txtProgressBar(max = n, style = 3) else NULL
}
