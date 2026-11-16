# ============================================================
# Kernel-smoothed local false discovery rate (local fdr) estimation,
# in the empirical-Bayes two-groups sense: assuming
# p_e ~ pi0 * f0(p) + pi1 * f1(p) with f0 uniform under the null,
# estimate the marginal density f(p) by kernel smoothing, estimate pi0
# by Storey's method, and set lfdr(p) = min(pi0 / f_hat(p), 1).
# ============================================================

#' Kernel-smooth a vector of (binned) counts
#'
#' @param counts Integer vector of bin counts.
#' @param K Smoothing window width (forced to odd).
#' @param kernel `"triangular"` (default) or `"uniform"`.
#'
#' @return Numeric vector, same length as `counts`.
#' @keywords internal
smooth_counts_discrete <- function(
  counts,
  K = 31,
  kernel = c("triangular", "uniform")
) {
  kernel <- match.arg(kernel)

  if (K %% 2 == 0) {
    K <- K + 1
  }

  n <- length(counts)
  h <- (K - 1) / 2
  out <- numeric(n)

  for (i in seq_len(n)) {
    lo <- max(1, i - h)
    hi <- min(n, i + h)
    d <- abs((lo:hi) - i)

    if (kernel == "uniform") {
      w <- rep(1, length(d))
    } else {
      w <- (h + 1) - d
    }

    w <- w / sum(w)
    out[i] <- sum(w * counts[lo:hi])
  }

  out
}

#' Local fdr from a Monte Carlo p-value on a discrete grid
#'
#' p-values from the simulation-based (pseudo-donor) null
#' (`evidence_p_method = "Q0_sim"`) live on the discrete grid
#' `{1/(B+1), ..., 1}`, so this estimates the p-value density on that
#' grid directly rather than treating it as continuous. Reached only via
#' `estimate_lfdr()`'s `evidence_p_method = "Q0_sim"` branch, which the
#' associated paper's analyses do not use (they use `"Q0_DL"`
#' throughout) -- kept as a legitimate, documented method option, but
#' not independently tested or validated by any analysis in that paper.
#'
#' @param p Numeric vector of p-values.
#' @param n_sim Number of null simulation replicates `B` (determines
#'   grid resolution as `n_sim + 1`); ignored if `M` is supplied.
#' @param M Explicit grid resolution, overriding `n_sim`.
#' @param K Smoothing window width, passed to [smooth_counts_discrete()].
#' @param kernel Smoothing kernel, passed to [smooth_counts_discrete()].
#' @param pi0_method `"storey"` (Storey's method, via `pi0_lambda`) or
#'   `"fixed"` (use `pi0_fixed` directly).
#' @param pi0_fixed Fixed `pi0` value, used when `pi0_method = "fixed"`.
#' @param pi0_lambda Storey's tuning parameter, used when
#'   `pi0_method = "storey"`.
#' @param isotonic `"pava"` (isotonic regression via `stats::isoreg()`),
#'   `"cummax"`, or `"none"`, enforcing `lfdr(p)` monotone in `p`.
#' @param eps Numerical floor/ceiling for clipping probabilities away
#'   from 0/1.
#'
#' @return A list with `p`, `p_raw`, `pi0`, `pi1`, `f_p`, `lfdr`,
#'   `rho_e` (= `1 - lfdr`), `discrete_grid` (grid-level diagnostics),
#'   and `settings`.
#' @keywords internal
estimate_discrete_lfdr <- function(
  p,
  n_sim = NULL,
  M = NULL,
  K = 31,
  kernel = c("triangular", "uniform"),
  pi0_method = c("storey", "fixed"),
  pi0_fixed = NULL,
  pi0_lambda = 0.5,
  isotonic = c("pava", "cummax", "none"),
  eps = 1e-12
) {

  kernel <- match.arg(kernel)
  pi0_method <- match.arg(pi0_method)
  isotonic <- match.arg(isotonic)

  p_raw <- p
  p <- as.numeric(p)
  ok <- is.finite(p)
  p <- pmin(pmax(p[ok], eps), 1 - eps)

  # Grid resolution: M overrides n_sim+1 when provided, decoupling density
  # estimation quality from simulation count.
  if (!is.null(M)) {
    n_grid <- as.integer(M)
  } else {
    if (is.null(n_sim)) {
      p_unique <- sort(unique(p))
      if (length(p_unique) < 2) stop("Either M or n_sim must be provided.")
      delta_hat <- min(diff(p_unique))
      n_sim <- round(1 / delta_hat - 1)
    }
    n_grid <- as.integer(n_sim) + 1
  }
  grid <- (1:n_grid) / n_grid
  delta <- 1 / n_grid
  bin_id <- round(p * n_grid)
  bin_id <- pmin(pmax(bin_id, 1), n_grid)
  counts <- tabulate(bin_id, nbins = n_grid)
  counts_smooth <- smooth_counts_discrete(
    counts = counts,
    K = K,
    kernel = kernel
  )

  prob_grid <- counts_smooth / sum(counts_smooth)
  f_grid <- pmax(prob_grid / delta, eps)
  if (pi0_method == "fixed") {
    if (is.null(pi0_fixed)) stop("pi0_fixed must be provided when pi0_method='fixed'.")
    pi0 <- pi0_fixed
  } else {
    pi0 <- mean(p > pi0_lambda, na.rm = TRUE) / (1 - pi0_lambda)
  }

  pi0 <- min(max(pi0, 0), 1)
  pi1 <- 1 - pi0
  lfdr_raw <- pmin(pi0 / pmax(f_grid, eps), 1)
  if (isotonic == "pava") {
    lfdr_grid <- stats::isoreg(grid, lfdr_raw)$yf
  } else if (isotonic == "cummax") {
    lfdr_grid <- cummax(lfdr_raw)
  } else {
    lfdr_grid <- lfdr_raw
  }

  lfdr_grid <- pmin(pmax(lfdr_grid, 0), 1)
  rho_grid <- 1 - lfdr_grid

  f_out <- rep(NA_real_, length(p_raw))
  lfdr_out <- rep(NA_real_, length(p_raw))
  rho_out <- rep(NA_real_, length(p_raw))

  f_out[ok] <- f_grid[bin_id]
  lfdr_out[ok] <- lfdr_grid[bin_id]
  rho_out[ok] <- rho_grid[bin_id]

  list(
    p = p,
    p_raw = p_raw,
    pi0 = pi0,
    pi1 = pi1,
    f_p = f_out,
    lfdr = lfdr_out,
    rho_e = rho_out,
    density = NULL,
    discrete_grid = data.frame(
      p_grid = grid,
      counts = counts,
      counts_smooth = counts_smooth,
      f_grid = f_grid,
      lfdr_raw = lfdr_raw,
      lfdr = lfdr_grid,
      rho = rho_grid
    ),
    settings = list(
      n_sim = n_sim,
      K = K,
      kernel = kernel,
      pi0_method = pi0_method,
      pi0_fixed = pi0_fixed,
      pi0_lambda = pi0_lambda,
      isotonic = isotonic,
      density_method = "discrete",
      eps = eps,
      m = length(p)
    )
  )

}

#' Local fdr from a continuous p-value
#'
#' Compatible with [estimate_discrete_lfdr()], but bins `p` onto a
#' fixed-resolution continuous grid instead of one sized to a known
#' number of null simulation replicates. Used for
#' `evidence_p_method = "Q0_DL"` (the asymptotic chi-squared p-value),
#' the setting used throughout the associated paper's analyses.
#'
#' @inheritParams estimate_discrete_lfdr
#' @param M Grid resolution (must be >= 10).
#'
#' @return Same structure as [estimate_discrete_lfdr()].
#' @keywords internal
estimate_continuous_lfdr <- function(
  p,
  M = 512,
  K = 31,
  kernel = c("triangular", "uniform"),
  pi0_method = c("storey", "fixed"),
  pi0_fixed = NULL,
  pi0_lambda = 0.5,
  isotonic = c("pava", "cummax", "none"),
  eps = 1e-12
) {

  kernel <- match.arg(kernel)
  pi0_method <- match.arg(pi0_method)
  isotonic <- match.arg(isotonic)

  p_raw <- p
  p <- as.numeric(p)
  ok <- is.finite(p)
  p <- pmin(pmax(p[ok], eps), 1 - eps)
  M <- as.integer(M)
  if (M < 10) stop("M must be >= 10.")
  delta <- 1 / M
  grid_mid <- ((1:M) - 0.5) / M
  bin_id <- ceiling(p * M)
  bin_id <- pmin(pmax(bin_id, 1), M)
  counts <- tabulate(bin_id, nbins = M)
  counts_smooth <- smooth_counts_discrete(
    counts = counts,
    K = K,
    kernel = kernel
  )
  prob_grid <- counts_smooth / sum(counts_smooth)
  f_grid <- pmax(prob_grid / delta, eps)
  if (pi0_method == "fixed") {
    if (is.null(pi0_fixed)) stop("pi0_fixed must be provided when pi0_method='fixed'.")
    pi0 <- pi0_fixed
  } else {
    pi0 <- mean(p > pi0_lambda, na.rm = TRUE) / (1 - pi0_lambda)
  }

  pi0 <- min(max(pi0, 0), 1)
  pi1 <- 1 - pi0
  f_p <- f_grid[bin_id]
  lfdr_raw_grid <- pmin(pi0 / pmax(f_grid, eps), 1)
  if (isotonic == "pava") {
    lfdr_grid <- stats::isoreg(
      grid_mid,
      lfdr_raw_grid
    )$yf
  } else if (isotonic == "cummax") {
    lfdr_grid <- cummax(lfdr_raw_grid)
  } else {
    lfdr_grid <- lfdr_raw_grid
  }

  lfdr_grid <- pmin(pmax(lfdr_grid, 0), 1)
  rho_grid <- 1 - lfdr_grid
  lfdr <- lfdr_grid[bin_id]
  rho_e <- rho_grid[bin_id]

  f_out <- rep(NA_real_, length(p_raw))
  lfdr_out <- rep(NA_real_, length(p_raw))
  rho_out <- rep(NA_real_, length(p_raw))

  f_out[ok] <- f_p
  lfdr_out[ok] <- lfdr
  rho_out[ok] <- rho_e

  list(
    p = p,
    p_raw = p_raw,
    pi0 = pi0,
    pi1 = pi1,
    f_p = f_out,
    lfdr = lfdr_out,
    rho_e = rho_out,
    density = NULL,
    discrete_grid = data.frame(
      p_grid = grid_mid,
      counts = counts,
      counts_smooth = counts_smooth,
      f_grid = f_grid,
      lfdr_raw = lfdr_raw_grid,
      lfdr = lfdr_grid,
      rho = rho_grid
    ),
    settings = list(
      M = M,
      K = K,
      kernel = kernel,
      pi0_method = pi0_method,
      pi0_fixed = pi0_fixed,
      pi0_lambda = pi0_lambda,
      isotonic = isotonic,
      density_method = "continuous",
      eps = eps,
      m = length(p)
    )
  )

}

#' Local fdr dispatcher
#'
#' Dispatches to [estimate_continuous_lfdr()] (for
#' `evidence_p_method = "Q0_DL"`) or [estimate_discrete_lfdr()] (for
#' `"Q0_sim"`).
#'
#' @inheritParams estimate_discrete_lfdr
#' @param evidence_p_method `"Q0_DL"` or `"Q0_sim"`, matching the
#'   p-value's origin (see [fit_EWVDS_decomposition()]).
#' @param M Grid resolution passed to whichever of
#'   [estimate_continuous_lfdr()] / [estimate_discrete_lfdr()] is used.
#'   Defaults to `NULL`: the `"Q0_DL"` branch then falls back to 512
#'   internally, while the `"Q0_sim"` branch lets `n_sim` determine the
#'   grid instead. If you pass `n_sim` directly to this function for the
#'   `"Q0_sim"` branch, leave `M` at its default (or explicitly `NULL`)
#'   -- a non-`NULL` `M` takes priority over `n_sim` in
#'   [estimate_discrete_lfdr()] and will silently override it.
#'
#' @return Same structure as [estimate_discrete_lfdr()].
#' @keywords internal
estimate_lfdr <- function(
  p,
  evidence_p_method = c("Q0_DL","Q0_sim"),
  n_sim = NULL,
  M = NULL,
  K = 31,
  kernel = c("triangular", "uniform"),
  pi0_method = c("storey", "fixed"),
  pi0_fixed = NULL,
  pi0_lambda = 0.5,
  isotonic = c("pava", "cummax", "none"),
  eps = 1e-12
) {
  evidence_p_method <- match.arg(evidence_p_method)
  if (evidence_p_method == "Q0_DL") {
    estimate_continuous_lfdr(
      p = p,
      M = if (is.null(M)) 512 else M,
      K = K,
      kernel = kernel,
      pi0_method = pi0_method,
      pi0_fixed = pi0_fixed,
      pi0_lambda = pi0_lambda,
      isotonic = isotonic,
      eps = eps
    )
  } else {
    estimate_discrete_lfdr(
      p = p,
      n_sim = n_sim,
      M = M,
      K = K,
      kernel = kernel,
      pi0_method = pi0_method,
      pi0_fixed = pi0_fixed,
      pi0_lambda = pi0_lambda,
      isotonic = isotonic,
      eps = eps
    )
  }
}
