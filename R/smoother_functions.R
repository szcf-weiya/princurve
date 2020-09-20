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
