# Brownian motion curves on the grid, one per column
BrownMat = function(N, refinement)
{
    mat = matrix(nrow = refinement, ncol = N)
    c = 1
    while(c <= N)
    {
        vec = sde::BM(N = refinement-1)
        mat[,c] = vec
        c = c+1
    }
    return(mat)
}

# Gaussian-type kernel for the FAR operators, scaled to unit operator norm
funKernel = function(ref)
{
  Mat = matrix(0, ref, ref)

  for(i in 1:ref)
  {
    for(j in 1:ref)
    {
      Mat[i,j] = exp(0.5 * ((i/ref)^2 + (j/ref)^2))
    }
  }

  A = Mat / ref

  op_norm = norm(A, type = "2")

  Mat = Mat / op_norm

  return(Mat)
}

# Apply an integral operator to curves on the grid
funIntegral = function(ref, Mat, X)
{
    Mat = Mat %*% X
    return(Mat/ref)
}

# FAR kernels for lags 1 to p
generate_AR_kernels = function(ref, p, ar_strength = 0.6,
                                lag_strength = c("equal", "decay"))
{
  lag_strength = match.arg(lag_strength)
  Phi_list = list()
  base_kernel = funKernel(ref)

  for(i in 1:p)
  {
    if(lag_strength == "decay") {
      Phi_list[[i]] = base_kernel * (ar_strength^i)
    } else {
      Phi_list[[i]] = base_kernel * (ar_strength / p)
    }
  }

  return(Phi_list)
}

# Run the FAR(p) recursion on the innovations
funARMat_p = function(refinement, eta_star_val, Phi_list)
{
  p = length(Phi_list)
  Tn = ncol(eta_star_val)

  res = matrix(0, refinement, Tn)

  res[,1:p] = eta_star_val[,1:p]

  for(t in (p+1):Tn)
  {
    ar_sum = rep(0, refinement)

    for(i in 1:p)
    {
      ar_sum = ar_sum + funIntegral(refinement, Phi_list[[i]], res[,(t-i)])
    }

    res[,t] = ar_sum + eta_star_val[,t]
  }

  return(res)
}

# Brownian-motion covariance kernel for the FMA operators, scaled to unit operator norm
MAphi = function(ref)
{
  Mat = matrix(0, ref, ref)

  for(i in 1:ref)
  {
    for(j in 1:ref)
    {
      Mat[i,j] = min(i,j) / ref
    }
  }

  A = Mat / ref

  norm = norm(A, type = "2")

  Mat = Mat / norm

  return(Mat)
}

# FMA kernels for lags 1 to q
generate_MA_kernels = function(ref, q, ma_strength = 0.6,
                                lag_strength = c("equal", "decay"))
{
  lag_strength = match.arg(lag_strength)
  Theta_list = list()
  base_kernel = MAphi(ref)

  for(j in 1:q)
  {
    if(lag_strength == "decay") {
      Theta_list[[j]] = base_kernel * (ma_strength^j)
    } else {
      Theta_list[[j]] = base_kernel * (ma_strength / q)
    }
  }

  return(Theta_list)
}

# Apply the FMA(q) filter to the innovations
funMAMat_q = function(refinement, eta_star_val, Theta_list)
{
  q = length(Theta_list)
  Tn = ncol(eta_star_val)

  res = matrix(0, refinement, Tn)

  for(t in (q+1):Tn)
  {
    ma_sum = eta_star_val[,t]

    for(j in 1:q)
    {
      ma_sum = ma_sum +
        funIntegral(refinement, Theta_list[[j]], eta_star_val[,(t-j)])
    }

    res[,t] = ma_sum
  }

  return(res)
}
