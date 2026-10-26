# Fit the FARMA operators and check their stability
operator_fit = function(Y, p, q, Npc) {
  est = farma_estimation(Y, p = p, q = q, m = NULL, Npc = Npc)
  fitted = est$fitted_original
  residuals = est$residuals
  dimnames(fitted) = dimnames(Y)
  dimnames(residuals) = dimnames(Y)
  stability = farma_fit_diagnostics(est, Y)
  stability$reasons = unstable_reasons(est, stability)
  structure(
    list(order = c(p = p, q = q), phi = est$phi, theta = est$theta, fitted = fitted, residuals = residuals),
    class = "farfima_operators",
    Npc = Npc,
    stability = stability,
    estimator = est
  )
}

# Operators as a list of matrices with their panel titles, AR lags first
operator_list = function(x) {
  p = if (is.null(x$phi)) 0 else dim(x$phi)[3]
  q = if (is.null(x$theta)) 0 else dim(x$theta)[3]
  list(
    p = p,
    q = q,
    matrices = c(lapply(seq_len(p), function(j) x$phi[, , j]), lapply(seq_len(q), function(j) x$theta[, , j])),
    titles = c(sprintf("AR operator, lag %d", seq_len(p)), sprintf("MA operator, lag %d", seq_len(q)))
  )
}

# Shift the columns of a matrix by k lags, padding with NA
lagMat = function(Z, k)
{
  if (k > 0) {
    do.call(what = cbind,
            args = c(rep(NA, times = k),
                     list(Z[,-(ncol(Z)-(k-1):0), drop=FALSE])))
  } else if (k < 0) {
    do.call(what = cbind,
            args = c(list(Z[,-(1:abs(k)), drop=FALSE]),
                     rep(NA, times = abs(k))))
  } else {
    Z
  }
}

# Shift the columns by k lags for the autoregressive regressors
lagMat_AR = function(Z, k)
{
  if (k > 0) {
    cbind(matrix(NA, nrow(Z), k),
          Z[, 1:(ncol(Z)-k), drop=FALSE])
  } else {
    Z
  }
}

# Functional Hannan-Rissanen estimation of the FARMA(p, q) operators, code by Thomas Kuenzer
farma_estimation = function(X, p = 1, q = 0, m = NULL, Npc = 1)
{
  n = ncol(X)
  d = nrow(X)

  if (p < 0 || q < 0) stop("p and q must be >= 0")
  if (p == 0 && q == 0) stop("At least one of p or q must be > 0")

  if (is.null(m)) m = min(n - 2, 15 + p + q)

  Xm = rowMeans(X)

  Xc = t(scale(t(X), center = TRUE, scale = FALSE))

  Xe = eigen(cov(t(Xc)))$vectors
  Npc = min(Npc, d)

  # Project the estimated operators onto the leading eigenfunctions
  reduce_score_operators = function(operator_score_array) {
    n_lags = dim(operator_score_array)[3]
    reduced = array(NA_real_, dim = c(Npc, Npc, n_lags))

    for (lag_i in seq_len(n_lags)) {
      score_operator = matrix(
        operator_score_array[, , lag_i],
        nrow = d,
        ncol = Npc
      )
      reduced[, , lag_i] =
        crossprod(Xe[, seq_len(Npc), drop = FALSE], score_operator)
    }

    reduced
  }

  # FAR(p): least squares on the lagged leading scores
  if (q == 0)
  {
    Xs = crossprod(Xe[, 1:Npc, drop = FALSE], Xc)
    Ws = do.call(rbind,
                  lapply(1:p, lagMat_AR,
                         Z = Xs[1:Npc, , drop = FALSE]))

    Y = Xc[, (p+1):n, drop = FALSE]
    Ws = Ws[, (p+1):n, drop = FALSE]

    D = Y %*% t(Ws) / (n - p)
    G = Ws %*% t(Ws) / (n - p)

    A_hat = D %*% solve(G)

    dim(A_hat) = c(d, Npc, p)

    phi_reduced = reduce_score_operators(A_hat)

    A_hat = apply(A_hat, 3,
                   function(M) M %*% t(Xe[, 1:Npc, drop = FALSE]))

    dim(A_hat) = c(d, d, p)

    fitted = matrix(NA_real_, d, n)
    resid = matrix(NA_real_, d, n)

    for (t in (p+1):n)
    {
      pred = 0
      for (j in 1:p) {
        pred = pred + A_hat[,,j] %*% Xc[, t-j]
      }

      fitted[, t] = pred
      resid[, t] = Xc[, t] - pred
    }

    fitted_original = fitted + Xm

    return(list(
      phi = A_hat,
      theta = NULL,
      phi_reduced = phi_reduced,
      theta_reduced = NULL,
      fitted = fitted,
      fitted_original = fitted_original,
      residuals = resid,
      mean = Xm
    ))
  }

  # Stage one: a long FAR(m) fit gives estimated innovations
  Xs = crossprod(Xe, Xc)

  Ws = do.call(rbind,
                lapply(1:m, lagMat, Z = Xs[1:Npc, , drop=FALSE]))

  D = cov(t(Xs), t(Ws), use = "pairwise.complete.obs")
  opshat = D %*% solve(cov(t(Ws), use = "pairwise.complete.obs"))

  epshats = Xs - opshat %*% Ws
  epshats[is.na(epshats)] = 0

  # Stage two: regress the curves on lagged innovations and lagged scores
  if (p == 0) {
    Ws = do.call(rbind,
                  lapply(1:q, lagMat, Z = epshats[1:Npc, , drop=FALSE]))
  } else {
    Ws = rbind(
      do.call(rbind, lapply(1:q, lagMat, Z = epshats[1:Npc, , drop=FALSE])),
      do.call(rbind, lapply(1:p, lagMat, Z = Xs[1:Npc, , drop=FALSE]))
    )
  }

  D = cov(t(Xc), t(Ws), use = "pairwise.complete.obs")
  opshat = D %*% solve(cov(t(Ws), use = "pairwise.complete.obs"))

  dim(opshat) = c(d, Npc, ncol(opshat) / Npc)
  opshat_reduced = reduce_score_operators(opshat)
  opshat = apply(opshat, 3, "%*%", t(Xe[, 1:Npc, drop=FALSE]))
  dim(opshat) = c(d, d, p + q)

  thetahat = opshat[, , 1:q, drop=FALSE]
  thetahat_reduced = opshat_reduced[, , 1:q, drop = FALSE]
  if (p > 0) {
    phihat = opshat[, , (q + 1):(q + p), drop=FALSE]
    phihat_reduced =
      opshat_reduced[, , (q + 1):(q + p), drop = FALSE]
  } else {
    phihat = NULL
    phihat_reduced = NULL
  }

  # In-sample one-step-ahead fit and residuals
  maxLag = max(p, q)
  eps_hat = matrix(0, nrow = d, ncol = n)
  fitted = matrix(NA_real_, nrow = d, ncol = n)

  for (t in (maxLag + 1):n) {
    ar_part = 0
    if (p > 0) {
      for (j in 1:p) {
        ar_part = ar_part + phihat[, , j] %*% Xc[, t - j]
      }
    }

    ma_part = 0
    if (q > 0) {
      for (k in 1:q) {
        ma_part = ma_part + thetahat[, , k] %*% eps_hat[, t - k]
      }
    }

    fitted[, t] = ar_part + ma_part

    eps_hat[, t] = Xc[, t] - fitted[, t]
  }

  eps_hat[, 1:maxLag] = NA
  fitted[, 1:maxLag] = NA

  fitted_original = fitted + Xm

  return(list(
    phi = phihat,
    theta = thetahat,
    phi_reduced = phihat_reduced,
    theta_reduced = thetahat_reduced,
    fitted = fitted,
    fitted_original = fitted_original,
    residuals = eps_hat,
    mean = Xm
  ))
}
