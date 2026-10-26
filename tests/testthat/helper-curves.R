# Simulated curves shared by the tests: 21 grid points by 160 curves
X = sim_farfima(n = 160, d = 0.3, p = 1, q = 0, no_grid = 21, seed = 1)

# Error condition of an unstable FARMA fit, used to mock failed forecast origins
unstable = function(radius) {
  structure(class = c("farma_unstable_fit", "error", "condition"),
            list(message = "unstable FARMA fit (MA part not invertible)", call = NULL))
}
