test_that("check_covariates passes silently when all covariates are present", {
  skip_if_not_installed("SingleCellExperiment")
  cd <- data.frame(cell_type = rep("T", 6), Donor = rep(c("d1", "d2"), 3),
                    log_umi = rnorm(6))
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = matrix(1:24, nrow = 4, ncol = 6)), colData = cd
  )
  expect_invisible(check_covariates(sce, c("cell_type", "Donor")))
})

make_tiny_marginal_fixture <- function(seed = 101, n_cells = 150) {
  set.seed(seed)
  log1p_tc <- rnorm(n_cells, mean = 7, sd = 0.3)
  mu1 <- exp(-2   + 0.8 * log1p_tc)
  mu2 <- exp(-1.5 + 0.6 * log1p_tc)
  counts <- rbind(
    gene1 = rnbinom(n_cells, mu = mu1, size = 6),
    gene2 = rnbinom(n_cells, mu = mu2, size = 8)
  )
  cells <- paste0("c", seq_len(n_cells))
  colnames(counts) <- cells
  cd <- data.frame(
    cell_type           = rep("T", n_cells),
    Donor               = factor(rep(c("d1", "d2", "d3"), length.out = n_cells)),
    log1p_total_counts  = log1p_tc,
    row.names           = cells
  )
  SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts), colData = cd
  )
}

test_that("sc3workflow_no_copula fits a real single gene and returns valid PIT scores", {
  skip_if_not_installed("scDesign3")
  skip_if_not_installed("SingleCellExperiment")
  sce   <- make_tiny_marginal_fixture()
  sce_g <- sce["gene1", ]

  res <- sc3workflow_no_copula(
    sce              = sce_g,
    celltype         = "cell_type",
    other_covariates = "log1p_total_counts",
    family_use       = "nb",
    mu_formula       = "log1p_total_counts",
    sigma_formula    = "1"
  )

  expect_true(res$model_info$engine %in% c("gam", "gamlss"))
  n <- ncol(sce_g)
  expect_length(res$Y, n)
  expect_length(res$mu_hat, n)
  expect_false(anyNA(res$mu_hat))
  expect_true(all(res$U_randomized >= 0 & res$U_randomized <= 1))
  expect_true(all(is.finite(res$Z_randomized)))
})

test_that("run_marginal_stage runs the real per-gene pipeline end-to-end and assembles output", {
  skip_if_not_installed("scDesign3")
  skip_if_not_installed("SingleCellExperiment")
  sce <- make_tiny_marginal_fixture()
  comb_path <- tempfile("marginal_stage_test_")
  on.exit(unlink(comb_path, recursive = TRUE), add = TRUE)

  res <- run_marginal_stage(
    sce           = sce,
    comb_path     = comb_path,
    gene_names    = rownames(sce),
    covariates    = "log1p_total_counts",
    mu_formula    = "log1p_total_counts",
    sigma_formula = "1",
    family_use    = "nb",
    donor_col     = "Donor",
    n_cores       = 1L,
    timeout_sec   = 120L
  )

  expect_true(file.exists(file.path(res$gene_path, "1_gene1.rds")))
  expect_true(file.exists(file.path(res$gene_path, "2_gene2.rds")))

  n_cells <- ncol(sce)
  for (v in c("U_randomized_raw", "U_randomized", "Z_randomized_raw",
              "Z_randomized", "mu_hat", "sigma_hat")) {
    csv_path <- file.path(comb_path, paste0(v, ".csv"))
    expect_true(file.exists(csv_path))
    M <- read.csv(csv_path, row.names = 1, check.names = FALSE)
    expect_equal(dim(M), c(n_cells, 2L))
    expect_setequal(colnames(M), c("gene1", "gene2"))
  }

  # donor_col is kept in the saved metadata even though it's not one of
  # `covariates`, so run_donor_split() can find it later -- see
  # run_marginal_stage()'s own roxygen docs for `donor_col`.
  meta <- readRDS(file.path(comb_path, "metadata.rds"))
  expect_true("Donor" %in% colnames(meta$cell_covariates))
})

test_that("run_marginal_stage's gene_names/covariates/mu_formula/sigma_formula defaults match an all-genes, log1p_total_counts-only fit", {
  skip_if_not_installed("scDesign3")
  skip_if_not_installed("SingleCellExperiment")
  sce <- make_tiny_marginal_fixture()
  comb_path <- tempfile("marginal_stage_defaults_test_")
  on.exit(unlink(comb_path, recursive = TRUE), add = TRUE)

  res <- run_marginal_stage(
    sce         = sce,
    comb_path   = comb_path,
    donor_col   = "Donor",
    n_cores     = 1L,
    timeout_sec = 120L
  )

  expect_true(file.exists(file.path(res$gene_path, "1_gene1.rds")))
  expect_true(file.exists(file.path(res$gene_path, "2_gene2.rds")))
  meta <- readRDS(file.path(comb_path, "metadata.rds"))
  expect_setequal(colnames(meta$cell_covariates)[colnames(meta$cell_covariates) != "Donor"],
                  "log1p_total_counts")
})

test_that("run_marginal_stage auto-creates the celltype column when sce doesn't have one", {
  skip_if_not_installed("scDesign3")
  skip_if_not_installed("SingleCellExperiment")
  sce <- make_tiny_marginal_fixture()
  SummarizedExperiment::colData(sce)$cell_type <- NULL  # no real celltype column
  comb_path <- tempfile("marginal_stage_celltype_autofill_")
  on.exit(unlink(comb_path, recursive = TRUE), add = TRUE)

  res <- run_marginal_stage(
    sce         = sce,
    comb_path   = comb_path,
    donor_col   = "Donor",
    n_cores     = 1L,
    timeout_sec = 120L
  )

  expect_true(file.exists(file.path(res$gene_path, "1_gene1.rds")))
  expect_true(file.exists(file.path(res$gene_path, "2_gene2.rds")))
})

test_that("run_marginal_stage errors instead of silently reusing a comb_path cached for a different sce", {
  skip_if_not_installed("scDesign3")
  skip_if_not_installed("SingleCellExperiment")
  comb_path <- tempfile("marginal_stage_mismatch_")
  on.exit(unlink(comb_path, recursive = TRUE), add = TRUE)

  sce1 <- make_tiny_marginal_fixture(seed = 201, n_cells = 60)
  run_marginal_stage(
    sce = sce1, comb_path = comb_path, gene_names = rownames(sce1),
    covariates = "log1p_total_counts", mu_formula = "log1p_total_counts",
    sigma_formula = "1", family_use = "nb",
    donor_col = "Donor", n_cores = 1L, timeout_sec = 120L
  )

  sce2 <- make_tiny_marginal_fixture(seed = 202, n_cells = 90)
  expect_error(
    run_marginal_stage(
      sce = sce2, comb_path = comb_path, gene_names = rownames(sce2),
      covariates = "log1p_total_counts", mu_formula = "log1p_total_counts",
      sigma_formula = "1", family_use = "nb",
      donor_col = "Donor", n_cores = 1L, timeout_sec = 120L
    ),
    "already has cached Stage 1 output"
  )
})

test_that("check_covariates errors with the missing covariate's name when a covariate is absent", {
  skip_if_not_installed("SingleCellExperiment")
  cd <- data.frame(cell_type = rep("T", 6), Donor = rep(c("d1", "d2"), 3))
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = matrix(1:24, nrow = 4, ncol = 6)), colData = cd
  )
  expect_error(check_covariates(sce, c("cell_type", "nonexistent_col")),
               "Missing covariates: nonexistent_col")
})

test_that("run_marginal_stage(seed =) is bit-reproducible within one session and across n_cores", {
  skip_if_not_installed("scDesign3")
  skip_if_not_installed("SingleCellExperiment")
  skip_on_os("windows")
  sce <- make_tiny_marginal_fixture()
  run_u <- function(seed, n_cores) {
    p <- tempfile("seedtest_")
    on.exit(unlink(p, recursive = TRUE), add = TRUE)
    suppressMessages(run_marginal_stage(
      sce, comb_path = p, gene_names = rownames(sce), donor_col = "Donor",
      covariates = "log1p_total_counts", mu_formula = "log1p_total_counts",
      celltype = "cell_type", n_cores = n_cores, timeout_sec = 300L,
      seed = seed
    ))
    as.matrix(utils::read.csv(file.path(p, "U_randomized_raw.csv"), row.names = 1))
  }
  a <- run_u(7, 1L)
  b <- run_u(7, 1L)
  c2 <- run_u(7, 2L)
  d <- run_u(8, 1L)
  expect_identical(a, b)
  expect_identical(a, c2)
  expect_false(identical(a, d))
  expect_error(run_marginal_stage(sce, comb_path = tempfile(), seed = 1,
                                  fix_rng_streams = FALSE),
               "requires fix_rng_streams")
})
