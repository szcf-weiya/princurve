#' Smoother functions
#'
#' Each of these functions have an interface \code{function(lambda, xj, ...)}, and
#' return smoothed values for xj. The output is expected to be ordered along an ordered lambda.
#' This means that the following is true:
#' \preformatted{
#' x <- runif(100)
#' y <- runif(100)
#' ord <- sample.int(100)
#' sfun <- smoother_functions[[1]]
#' all(sfun(x, y) == sfun(x[ord], y[ord]))
#' }
#'
tophist = function(x, binwidth = 0.02) {
  # suppose has been sorted
  # x = sort(x)
  k = 1
  left_bin = min(x)
  counts = 0
  for (i in 1:length(x)) {
    if (left_bin[k] <= x[i] && x[i] < left_bin[k] + binwidth)
      counts[k] = counts[k] + 1
    else {
      tmp_left_bin = left_bin[k] + binwidth
      while (TRUE) {
        if (tmp_left_bin <= x[i] && x[i] < tmp_left_bin + binwidth) {
          k = k + 1
          left_bin = c(left_bin, tmp_left_bin)
          counts = c(counts, 1)
          break
        } else {
          tmp_left_bin = tmp_left_bin + binwidth
        }
      }
    }
  }
  left_bin[which.max(counts)] + binwidth / 2
}
sol_cubic_monotone_smooth_spline = function(y, B, Omega, alpha = 1, lambda = 1.0) {
  J = ncol(B)
  x = CVXR::Variable(J)
  A = diag(J) * alpha
  diag(A[1:J-1, 2:J]) = -1 * alpha
  b = numeric(J-1)
  Omega = (Omega + t(Omega)) / 2
  eig = eigen(Omega)
  d = eig$values
  if (sum(d < 0)) {
    cat("negative eigenvalues: ", d[d < 0], "\n")
    d[d < 0] = 0
    #Omega = eig$vectors %*% diag(d) %*% t(eig$vectors)
  }
  V = t(eig$vectors %*% diag(sqrt(d))) # Omega = PDP' = PD^{1/2} D^{1/2}P' := V'V
#  obj = CVXR::Minimize(CVXR::sum_squares(B %*% x - y) + lambda * CVXR::quad_form(x, Omega))
  obj = CVXR::Minimize(CVXR::sum_squares(B %*% x - y) + lambda * CVXR::sum_squares(V %*% x))
#  obj = CVXR::Minimize(CVXR::sum_squares(B %*% x - y))
  # alpha = 1: increasing
  # alpha = -1: decreasing
#  conds = list(A[1:J-1,] %*% x <= b,  CVXR::quad_form(x, Omega)<= 100)
  conds = list(A[1:J-1,] %*% x <= b)
  sol = CVXR::solve(CVXR::Problem(obj, conds))
  if (sol$status == "optimal")
    return(sol$getValue(x))
  else
    return(NULL)
}
smoother_functions <- list(
  smooth_spline = function(lambda, xj, ..., df = 5) {
    ord <- order(lambda)
    lambda <- lambda[ord]
    xj <- xj[ord]
    fit <- stats::smooth.spline(lambda, xj, ..., df = df, keep.data = FALSE)
    # cat(fit$df, '\n')
    stats::predict(fit, x = lambda)$y
  },

  # doc of smooth_spline said
  # "If spar and lambda are missing or NULL,
  #  the value of df is used to determine the degree of smoothing.
  #  If df is missing as well, leave-one-out cross-validation
  #  (ordinary or ‘generalized’ as determined by cv) is used to determine λ (penalty parameter)."
  cv_smooth_spline = function(lambda, xj, ...) {
    ord <- order(lambda)
    lambda <- lambda[ord]
    xj <- xj[ord]
    fit <- stats::smooth.spline(lambda, xj, ..., keep.data = FALSE)
    # cat(fit$df, '\n')
    stats::predict(fit, x = lambda)$y
  },

  monotone_cubic_spline = function(lambda, xj, ..., J  = 10) {
    ord <- order(lambda)
    lambda <- lambda[ord]
    xj <- xj[ord]
    H = splines::bs(lambda, df = J, intercept = TRUE)
    if (sum(xj > -0.1) > 60) { # dop
      A = diag(J)
      diag(A[1:J-1, 2:J]) = -1
      b = numeric(J-1)
      beta = lsei::lsi(H, xj, e = A[1:J-1,], f = b)
      predict(H, x = lambda) %*% beta
    } else {
      beta = lm(xj ~ 0 + H)$coefficients
      predict(H, x = lambda) %*% beta
    }
  },

  double_monotone_cubic_spline = function(lambda, xj, ..., J  = 10) {
    ord <- order(lambda)
    lambda <- lambda[ord]
    xj <- xj[ord]
    H = splines::bs(lambda, df = J, intercept = TRUE)
    A = diag(J)
    diag(A[1:J-1, 2:J]) = -1
    b = numeric(J-1)
    if (sum(xj > -0.1) > 60) { # dop
      beta = lsei::lsi(H, xj, e = A[1:J-1,], f = b)
    } else {
      beta = lsei::lsi(H, xj, e = -A[1:J-1,], f = b)
    }
    predict(H, x = lambda) %*% beta
  },

  double_monotone_cubic_smooth_spline = function(lambda, xj, ..., mu = 1.0, by=5) {
    ord <- order(lambda)
    lambda <- lambda[ord]
    xj <- xj[ord]
    n = length(lambda)
    bbasis = fda::create.bspline.basis(breaks = unique(lambda[seq(1, n, by=by)]), norder = 4)
    J = bbasis$nbasis
    Omega = fda::eval.penalty(bbasis, 2)
    B = fda::eval.basis(lambda, bbasis)
    if (sum(xj > -0.1) > 60) #dop
      alpha = -1
    else
      alpha = 1
    gamma = sol_cubic_monotone_smooth_spline(xj, B, Omega, alpha = alpha, lambda = mu)
    if (is.null(gamma)){
      cat("No optimal solution, and use double monotone cubic spline")
      double_monotone_cubic_spline(lambda, xj, ...)
    } else {
      B %*% gamma
    }
  },

  w_smooth_spline = function(lambda, xj, ...) {
    ord <- order(lambda)
    lambda <- lambda[ord]
    xj <- xj[ord]
    lambda0 = tophist(lambda)
    fit <- stats::smooth.spline(lambda, xj, 1 / (abs(lambda - lambda0) + 1), ..., keep.data = FALSE)
    stats::predict(fit, x = lambda)$y
  },

  log_cv_smooth_spline = function(lambda, xj, ...) {
    ord <- order(lambda)
    lambda <- lambda[ord]
    xj <- xj[ord]
    fit <- stats::smooth.spline(log(lambda - lambda[1] + 1.1), xj, ..., keep.data = FALSE)
    # cat(fit$df, '\n')
    stats::predict(fit, x = log(lambda - lambda[1] + 1.1))$y
  },


  lowess = function(lambda, xj, ...) {
    stats::lowess(lambda, xj, ...)$y
  },

  periodic_lowess = function(lambda, xj, ...) {
    periodic_lowess(lambda, xj, ...)$y
  }
)
