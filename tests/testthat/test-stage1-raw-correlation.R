make_tiny_donor_split_fixture <- function(seed = 21) {
  set.seed(seed)
  comb_path <- tempfile("donor_split_test_")
  dir.create(comb_path)

  genes <- c("g1", "g2", "g3")
  cells <- paste0("c", 1:12)
  donor_vec <- rep(c("d1", "d2"), each = 6)
  names(donor_vec) <- cells

  U <- matrix(runif(12 * 3), nrow = 12, ncol = 3, dimnames = list(cells, genes))
  Z <- qnorm(U)
  write.csv(U, file.path(comb_path, "U_randomized.csv"), row.names = TRUE)
  write.csv(Z, file.path(comb_path, "Z_randomized.csv"), row.names = TRUE)
  saveRDS(list(cell_covariates = data.frame(Donor = donor_vec, row.names = cells)),
          file.path(comb_path, "metadata.rds"))

  list(comb_path = comb_path, genes = genes, cells = cells, donor_vec = donor_vec, Z = Z)
}

test_that("run_donor_split writes Sigma_global.rds with diag 0 and the right dimensions", {
  fx <- make_tiny_donor_split_fixture()
  on.exit(unlink(fx$comb_path, recursive = TRUE), add = TRUE)

  res <- run_donor_split(comb_path = fx$comb_path, donor_col = "Donor", n_cores = 1)

  expect_setequal(res$donors, c("d1", "d2"))
  sig_path <- file.path(fx$comb_path, "Sigma_global.rds")
  expect_true(file.exists(sig_path))
  Sigma_global <- readRDS(sig_path)
  expect_equal(dim(Sigma_global), c(length(fx$genes), length(fx$genes)))
  # diag is explicitly zeroed out here (edge-level global correlation is
  # read off the off-diagonal only; unlike a donor's own Sigma_Z.csv,
  # which keeps the usual diag = 1 from a plain cor() call).
  expect_equal(unname(diag(Sigma_global)), rep(0, length(fx$genes)))
})

test_that("run_donor_split's per-donor Sigma_Z.csv matches a direct cor() on that donor's own cells", {
  fx <- make_tiny_donor_split_fixture()
  on.exit(unlink(fx$comb_path, recursive = TRUE), add = TRUE)

  run_donor_split(comb_path = fx$comb_path, donor_col = "Donor", n_cores = 1)

  for (d in c("d1", "d2")) {
    donor_dir <- file.path(fx$comb_path, "donors", d)
    expect_true(file.exists(file.path(donor_dir, "U.csv")))
    expect_true(file.exists(file.path(donor_dir, "Z.csv")))
    expect_true(file.exists(file.path(donor_dir, "Sigma_Z.csv")))

    Sigma_d <- as.matrix(read.csv(file.path(donor_dir, "Sigma_Z.csv"), row.names = 1, check.names = FALSE))
    Z_d_direct <- fx$Z[fx$donor_vec == d, , drop = FALSE]
    expect_equal(unname(Sigma_d), unname(cor(Z_d_direct)))
  }
})

test_that("run_donor_split skips recomputing an existing Sigma_global.rds", {
  fx <- make_tiny_donor_split_fixture()
  on.exit(unlink(fx$comb_path, recursive = TRUE), add = TRUE)

  run_donor_split(comb_path = fx$comb_path, donor_col = "Donor", n_cores = 1)
  sig_path <- file.path(fx$comb_path, "Sigma_global.rds")
  first_mtime <- file.info(sig_path)$mtime

  Sys.sleep(1.1)  # ensure a re-write, if it happened, would have a detectably later mtime
  run_donor_split(comb_path = fx$comb_path, donor_col = "Donor", n_cores = 1)
  second_mtime <- file.info(sig_path)$mtime

  expect_equal(first_mtime, second_mtime)
})

test_that("run_donor_split's Sigma_global.rds is cor() of the same Z used for donors", {
  fx <- make_tiny_donor_split_fixture()
  on.exit(unlink(fx$comb_path, recursive = TRUE), add = TRUE)

  run_donor_split(comb_path = fx$comb_path, donor_col = "Donor", n_cores = 1)

  Sigma_global <- readRDS(file.path(fx$comb_path, "Sigma_global.rds"))
  expected <- cor(fx$Z)
  diag(expected) <- 0
  expect_equal(unname(Sigma_global), unname(expected))
})

test_that("run_donor_split and run_coexp_pipeline warn on the deprecated clip arguments", {
  fx <- make_tiny_donor_split_fixture()
  on.exit(unlink(fx$comb_path, recursive = TRUE), add = TRUE)

  expect_warning(
    run_donor_split(comb_path = fx$comb_path, donor_col = "Donor", eps = 1e-4, n_cores = 1),
    "deprecated"
  )
  expect_warning(
    tryCatch(run_coexp_pipeline(sce = NULL, comb_path = tempfile(), donor_col = "Donor",
                                eps_pit = 1e-4),
             error = function(e) NULL),
    "deprecated"
  )
})
