test_that("smooth_counts_discrete leaves a constant vector unchanged", {
  expect_equal(smooth_counts_discrete(rep(5, 10), K = 5, kernel = "uniform"), rep(5, 10))
  expect_equal(smooth_counts_discrete(rep(5, 10), K = 5, kernel = "triangular"), rep(5, 10))
})

test_that("smooth_counts_discrete matches an independent re-derivation of the weighted average", {
  set.seed(1)
  counts <- sample(0:5, 20, replace = TRUE)
  h <- 2  # K = 5 -> h = (K-1)/2

  manual_smooth <- function(counts, h, kernel) {
    n <- length(counts)
    out <- numeric(n)
    for (i in seq_len(n)) {
      lo <- max(1, i - h); hi <- min(n, i + h)
      d <- abs((lo:hi) - i)
      w <- if (kernel == "uniform") rep(1, length(d)) else (h + 1) - d
      w <- w / sum(w)
      out[i] <- sum(w * counts[lo:hi])
    }
    out
  }

  expect_equal(smooth_counts_discrete(counts, K = 5, kernel = "uniform"),
               manual_smooth(counts, h, "uniform"))
  expect_equal(smooth_counts_discrete(counts, K = 5, kernel = "triangular"),
               manual_smooth(counts, h, "triangular"))

  # boundary point (window truncated, not reflected/padded)
  expect_equal(smooth_counts_discrete(counts, K = 5, kernel = "triangular")[1],
               manual_smooth(counts, h, "triangular")[1])
})

test_that("estimate_continuous_lfdr: an all-null p-value vector gives pi0 near 1 and lfdr near 1 everywhere", {
  set.seed(123)
  p_null <- runif(3000)
  fit <- estimate_continuous_lfdr(p_null, M = 200, pi0_lambda = 0.5)
  expect_gt(fit$pi0, 0.95)
  expect_gt(min(fit$lfdr), 0.9)
})

test_that("estimate_continuous_lfdr: an all-signal p-value vector gives pi0 near 0", {
  set.seed(123)
  p_signal <- rbeta(3000, 0.3, 8)  # concentrated near 0
  fit <- estimate_continuous_lfdr(p_signal, M = 200, pi0_lambda = 0.5)
  expect_lt(fit$pi0, 0.05)
})

test_that("estimate_continuous_lfdr: lfdr is monotone non-decreasing in p under the default isotonic regression", {
  set.seed(123)
  p <- rbeta(3000, 0.3, 8)
  fit <- estimate_continuous_lfdr(p, M = 200, pi0_lambda = 0.5)
  o <- order(p)
  expect_true(all(diff(fit$lfdr[o]) >= -1e-9))
})

test_that("estimate_continuous_lfdr: pi0_method = 'fixed' passes pi0_fixed straight through", {
  set.seed(123)
  p <- runif(1000)
  fit <- estimate_continuous_lfdr(p, M = 200, pi0_method = "fixed", pi0_fixed = 0.7)
  expect_equal(fit$pi0, 0.7)
})

test_that("estimate_discrete_lfdr: an all-null discrete-grid p-value vector gives pi0 near 1 and lfdr near 1", {
  set.seed(99)
  B <- 999
  p_null_disc <- (1 + sample(0:B, 3000, replace = TRUE)) / (B + 1)
  fit <- estimate_discrete_lfdr(p_null_disc, n_sim = B, pi0_lambda = 0.5)
  expect_gt(fit$pi0, 0.9)
  expect_gt(min(fit$lfdr), 0.85)
})

test_that("estimate_discrete_lfdr: M, when supplied, overrides n_sim for grid resolution", {
  set.seed(99)
  B <- 999
  p_null_disc <- (1 + sample(0:B, 1000, replace = TRUE)) / (B + 1)
  fit_n_sim_only <- estimate_discrete_lfdr(p_null_disc, n_sim = B)
  fit_M_override <- estimate_discrete_lfdr(p_null_disc, n_sim = B, M = 50)
  expect_equal(nrow(fit_n_sim_only$discrete_grid), B + 1)
  expect_equal(nrow(fit_M_override$discrete_grid), 50)
})

test_that("estimate_lfdr dispatches Q0_DL to estimate_continuous_lfdr with the same default M = 512 fallback", {
  set.seed(1)
  p <- runif(500)
  expect_identical(estimate_lfdr(p, evidence_p_method = "Q0_DL"),
                    estimate_continuous_lfdr(p, M = 512))
})

test_that("estimate_lfdr dispatches Q0_sim to estimate_discrete_lfdr, using n_sim when M is left NULL", {
  set.seed(1)
  p <- (1 + sample(0:99, 500, replace = TRUE)) / 100
  expect_identical(estimate_lfdr(p, evidence_p_method = "Q0_sim", n_sim = 99),
                    estimate_discrete_lfdr(p, n_sim = 99))
})
