test_that("clip_cor clips values outside (-1+eps, 1-eps)", {
  expect_equal(clip_cor(c(-2, -1, 0, 1, 2), eps = 1e-6),
               c(-1 + 1e-6, -1 + 1e-6, 0, 1 - 1e-6, 1 - 1e-6))
  expect_equal(clip_cor(0.5), 0.5)
})

test_that("fisher_z matches atanh on clipped input", {
  expect_equal(fisher_z(0.5), atanh(0.5))
  expect_equal(fisher_z(0), 0)
  # saturated correlations are clipped before atanh, not passed straight through
  expect_equal(fisher_z(1, eps = 1e-6), atanh(1 - 1e-6))
  expect_equal(fisher_z(-1, eps = 1e-6), atanh(-1 + 1e-6))
})

test_that("inv_fisher_z is tanh and round-trips with fisher_z", {
  expect_equal(inv_fisher_z(0), 0)
  r <- c(-0.9, -0.3, 0, 0.3, 0.9)
  expect_equal(inv_fisher_z(fisher_z(r)), r, tolerance = 1e-8)
})

test_that("get_upper_idx returns exactly the strict upper triangle, in column-major order", {
  Sigma <- matrix(c(1, 0.2, 0.3,
                     0.2, 1, 0.4,
                     0.3, 0.4, 1), nrow = 3)
  idx <- get_upper_idx(Sigma)
  expect_equal(dim(idx), c(3L, 2L))
  expect_equal(unname(idx), matrix(c(1, 1, 2,
                                      2, 3, 3), nrow = 3))
  expect_equal(Sigma[cbind(idx[, 1], idx[, 2])], c(0.2, 0.3, 0.4))
})

test_that("edge_to_symmetric_mat round-trips with get_upper_idx", {
  Sigma <- matrix(c(1, 0.2, 0.3,
                     0.2, 1, 0.4,
                     0.3, 0.4, 1), nrow = 3)
  idx <- get_upper_idx(Sigma)
  edge_values <- Sigma[cbind(idx[, 1], idx[, 2])]
  rebuilt <- edge_to_symmetric_mat(edge_values, idx, diag_value = 1)
  expect_equal(rebuilt, Sigma)
})

test_that("edge_to_symmetric_mat sets the requested diagonal value and gene names", {
  idx <- matrix(c(1, 1, 2, 2, 3, 3), nrow = 3)
  rebuilt <- edge_to_symmetric_mat(c(0.1, 0.2, 0.3), idx,
                                    gene_names = c("A", "B", "C"),
                                    diag_value = 0)
  expect_equal(unname(diag(rebuilt)), c(0, 0, 0))
  expect_equal(rownames(rebuilt), c("A", "B", "C"))
  expect_equal(colnames(rebuilt), c("A", "B", "C"))
  expect_true(isSymmetric(rebuilt))
})

# ---- ensure_psd_correlation ----

# A symmetric, unit-diagonal matrix with entries in (-1, 1) that is not PSD:
# genes 1-2 and 1-3 strongly positive but 2-3 strongly negative.
make_non_psd_cor <- function() {
  S <- matrix(c( 1,   0.9,  0.9,
                 0.9, 1,   -0.9,
                 0.9, -0.9, 1), nrow = 3)
  dimnames(S) <- list(paste0("g", 1:3), paste0("g", 1:3))
  S
}

test_that("ensure_psd_correlation returns an already-PSD matrix unchanged", {
  S <- matrix(c(1, 0.2, 0.3,
                0.2, 1, 0.4,
                0.3, 0.4, 1), nrow = 3)
  res <- ensure_psd_correlation(S)
  expect_false(res$projected)
  expect_identical(res$mat, S)
  expect_equal(res$max_abs_change, 0)
  expect_equal(res$min_eig_before, res$min_eig_after)
})

test_that("ensure_psd_correlation projects a non-PSD matrix to a valid correlation matrix", {
  S <- make_non_psd_cor()
  expect_lt(min(eigen(S, only.values = TRUE)$values), 0)

  res <- ensure_psd_correlation(S)
  expect_true(res$projected)
  expect_lt(res$min_eig_before, 0)
  expect_gte(res$min_eig_after, -1e-8)
  expect_equal(unname(diag(res$mat)), rep(1, 3))
  expect_equal(res$mat, t(res$mat))
  expect_identical(dimnames(res$mat), dimnames(S))
  expect_true(is.matrix(res$mat))  # base matrix, not a Matrix-package S4 object
  expect_gt(res$max_abs_change, 0)
})

test_that("ensure_psd_correlation treats eigenvalues within -tol as PSD", {
  S <- make_non_psd_cor()
  ev <- min(eigen(S, only.values = TRUE)$values)
  expect_false(ensure_psd_correlation(S, tol = abs(ev) + 1e-6)$projected)
})
