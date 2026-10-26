#' Simulate a FARFIMA(p, d, q) process
#'
#' Simulates `n` curves from a FARFIMA(p, d, q) model on `no_grid` equally
#' spaced points of \[0, 1\]. The innovations are Brownian motion curves. They
#' are fractionally integrated with memory parameter `d`, passed through the
#' MA operators and then through the AR recursion. The first 100 curves are
#' dropped as burn-in.
#'
#' The AR kernel is \eqn{\exp\{(u^2 + v^2)/2\}} and the MA kernel is
#' \eqn{\min(u, v)}. Each is scaled to operator norm one and multiplied by
#' `ar_strength` or `ma_strength`.
#'
#' @param n Number of curves.
#' @param d Memory parameter.
#' @param p,q AR and MA orders, not both zero.
#' @param no_grid Number of grid points.
#' @param seed Optional seed.
#' @param ar_strength,ma_strength Operator strengths \eqn{\psi} and \eqn{\nu}.
#' @param lag_strength `"equal"` splits the strength equally across lags, and
#'   `"decay"` lets it decay geometrically.
#' @return A matrix of curves, `no_grid` by `n`. The innovations are kept in
#'   `attr(, "innovations")`.
#' @seealso [split_curves()], [farfima()]
#' @examples
#' sim = sim_farfima(n = 200, d = 0.2, p = 1, q = 1, seed = 1)
#' dim(sim)
#' @export
sim_farfima = function(
    n = 500,
    d = 0.1,
    p = 1,
    q = 0,
    no_grid = 101,
    seed = NULL,
    ar_strength = 0.6,
    ma_strength = 0.6,
    lag_strength = c("equal", "decay")
)
{
  lag_strength = match.arg(lag_strength)

  if (!is.null(seed)) {
    withr::local_seed(seed)
  }

  # Fractionally integrated Brownian motion innovations, with 100 extra curves as burn-in
  eps = BrownMat(N = n + 100, refinement = no_grid)

  eta_star = t(FCVAR::FracDiff(
    x = t(eps),
    d = -d
  ))

  # Pass the innovations through the FAR, FMA or FARMA operators
  if(p != 0 & q == 0)
  {
    Phi_list = generate_AR_kernels(no_grid, p, ar_strength, lag_strength)
    sim_X = funARMat_p(no_grid, eta_star, Phi_list)

  } else if(p == 0 & q != 0)
  {
    Theta_list = generate_MA_kernels(no_grid, q, ma_strength, lag_strength)
    sim_X = funMAMat_q(no_grid, eta_star, Theta_list)

  } else if(p != 0 & q != 0)
  {
    Theta_list = generate_MA_kernels(no_grid, q, ma_strength, lag_strength)
    Phi_list = generate_AR_kernels(no_grid, p, ar_strength, lag_strength)

    eta_star_MA = funMAMat_q(no_grid, eta_star, Theta_list)
    sim_X = funARMat_p(no_grid, eta_star_MA, Phi_list)
  }
  else
  {
    stop("`p` and `q` cannot both be zero.", call. = FALSE)
  }

  sim_X_record = sim_X[, (ncol(sim_X) - (n - 1)):ncol(sim_X)]
  eta_record = eta_star[, (ncol(eta_star) - (n - 1)):ncol(eta_star)]

  structure(sim_X_record, innovations = eta_record)
}
