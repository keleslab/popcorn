test_that("generate_pseudodonor_indices: resampling mode draws the right length per donor, within range", {
  set.seed(1)
  dn <- c(d1 = 5, d2 = 8, d3 = 3)
  idx_list <- generate_pseudodonor_indices(dn, n_cell_total = 20, mode = "resampling", replace = TRUE)
  expect_equal(names(idx_list), names(dn))
  expect_equal(lengths(idx_list), unname(dn), ignore_attr = "names")
  expect_true(all(unlist(idx_list) >= 1 & unlist(idx_list) <= 20))
})

test_that("generate_pseudodonor_indices: permutation mode gives non-overlapping index sets", {
  set.seed(1)
  dn <- c(d1 = 5, d2 = 8, d3 = 3)
  idx_list <- generate_pseudodonor_indices(dn, n_cell_total = 20, mode = "permutation")
  expect_equal(lengths(idx_list), unname(dn), ignore_attr = "names")
  all_idx <- unlist(idx_list)
  expect_equal(length(all_idx), length(unique(all_idx)))
})

test_that("generate_pseudodonor_indices: permutation mode errors when donors need more cells than available", {
  expect_error(
    generate_pseudodonor_indices(c(d1 = 15, d2 = 15), n_cell_total = 20, mode = "permutation"),
    "sum\\(donor_ncells\\)"
  )
})

test_that("compute_pseudodonor_Zmat: correct dimensions and Y_null_mat = Z_null_mat - Z_g_edge exactly", {
  set.seed(2)
  Z_cell_mat <- matrix(rnorm(20 * 4), nrow = 20, ncol = 4,
                        dimnames = list(NULL, c("g1", "g2", "g3", "g4")))
  idx <- get_upper_idx(matrix(0, 4, 4))
  Z_g_edge <- fisher_z(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6))
  dn <- c(d1 = 5, d2 = 8, d3 = 3)
  pseudo_idx_list <- generate_pseudodonor_indices(dn, n_cell_total = 20, mode = "resampling", replace = FALSE)

  res <- compute_pseudodonor_Zmat(Z_cell_mat, Z_g_edge, pseudo_idx_list, idx)

  expect_equal(dim(res$Z_null_mat), c(nrow(idx), length(dn)))
  expect_equal(colnames(res$Z_null_mat), names(dn))
  expect_equal(res$Y_null_mat, sweep(res$Z_null_mat, 1, Z_g_edge, "-"))
})

test_that("simulate_Q0_null_pseudodonor: correct shape, finite and non-negative, reproducible under a fixed seed", {
  set.seed(3)
  Z_cell_mat <- matrix(rnorm(30 * 4), nrow = 30, ncol = 4,
                        dimnames = list(NULL, c("g1", "g2", "g3", "g4")))
  idx <- get_upper_idx(matrix(0, 4, 4))
  Z_g_edge <- fisher_z(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6))
  dn <- c(d1 = 10, d2 = 10, d3 = 10)
  W0 <- matrix(rep(dn - 3, each = nrow(idx)), nrow = nrow(idx), ncol = length(dn),
               dimnames = list(NULL, names(dn)))

  res <- simulate_Q0_null_pseudodonor(Z_cell_mat, Z_g_edge, dn, idx, W0,
                                       n_sim = 5, mode = "resampling",
                                       seed = 10, n_cores_sim = 1)
  expect_equal(dim(res$Q0_null_mat), c(nrow(idx), 5L))
  expect_true(all(is.finite(res$Q0_null_mat)))
  expect_true(all(res$Q0_null_mat >= 0))

  res2 <- simulate_Q0_null_pseudodonor(Z_cell_mat, Z_g_edge, dn, idx, W0,
                                        n_sim = 5, mode = "resampling",
                                        seed = 10, n_cores_sim = 1)
  expect_identical(res$Q0_null_mat, res2$Q0_null_mat)
})

test_that("simulate_Q0_null_pseudodonor is reproducible under a fixed seed on several cores", {
  skip_on_os("windows")
  set.seed(3)
  Z_cell_mat <- matrix(rnorm(30 * 4), nrow = 30, ncol = 4,
                        dimnames = list(NULL, c("g1", "g2", "g3", "g4")))
  idx <- get_upper_idx(matrix(0, 4, 4))
  Z_g_edge <- fisher_z(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6))
  dn <- c(d1 = 10, d2 = 10, d3 = 10)
  W0 <- matrix(rep(dn - 3, each = nrow(idx)), nrow = nrow(idx), ncol = length(dn),
               dimnames = list(NULL, names(dn)))
  kind_before <- RNGkind()[1]
  f <- function() simulate_Q0_null_pseudodonor(Z_cell_mat, Z_g_edge, dn, idx, W0,
                                              n_sim = 8, mode = "resampling",
                                              seed = 10, n_cores_sim = 2)$Q0_null_mat
  expect_identical(f(), f())
  expect_identical(RNGkind()[1], kind_before)
})

test_that("pseudo-donor null defaults to disjoint (permutation) mode everywhere", {
  dn <- c(d1 = 5, d2 = 7, d3 = 3)
  set.seed(1)
  idx_list <- generate_pseudodonor_indices(dn, n_cell_total = 20)
  expect_equal(anyDuplicated(unlist(idx_list)), 0L)
  expect_equal(eval(formals(run_coexp_pipeline)$pseudo_donor_mode), "permutation")
  expect_equal(eval(formals(run_ewvds_stage)$pseudo_donor_mode), "permutation")
  expect_equal(eval(formals(run_EWVDS_pipeline)$pseudo_donor_mode)[1], "permutation")
  expect_equal(eval(formals(fit_EWVDS_decomposition)$pseudo_donor_mode)[1], "permutation")
  expect_equal(eval(formals(simulate_Q0_null_pseudodonor)$mode)[1], "permutation")
})

test_that("internal Stage 2 defaults match the exported entry points", {
  for (f in list(run_EWVDS_pipeline, fit_EWVDS_decomposition)) {
    fm <- formals(f)
    expect_false(eval(fm$compute_Q0_sim))
    expect_equal(eval(fm$null_method)[1], "pseudo_donor")
    expect_equal(eval(fm$excess_variance_method)[1], "DL_tau")
    expect_equal(eval(fm$evidence_p_method)[1], "Q0_DL")
  }
  expect_equal(eval(formals(run_coexp_pipeline)$pipeline_tag),
               eval(formals(run_ewvds_stage)$pipeline_tag))
})
