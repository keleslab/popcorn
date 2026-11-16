test_that("bootstrap_edge_variance: correct length, non-negative, reproducible under a fixed seed", {
  set.seed(5)
  Z_donor <- matrix(rnorm(50 * 4), nrow = 50, ncol = 4)
  idx <- get_upper_idx(matrix(0, 4, 4))

  v1 <- bootstrap_edge_variance(Z_donor, idx, n_boot = 30, seed = 7)
  expect_equal(length(v1), nrow(idx))
  expect_true(all(v1 >= 0))

  v2 <- bootstrap_edge_variance(Z_donor, idx, n_boot = 30, seed = 7)
  expect_identical(v1, v2)
})

test_that("remove_EWVDS_pcs exactly removes an exact rank-1 deviation matrix", {
  u <- c(1, 2, -1, 0.5, 3)
  v <- c(2, -1, 0.5)
  Y_rank1 <- outer(u, v)
  colnames(Y_rank1) <- c("d1", "d2", "d3")

  resid <- remove_EWVDS_pcs(Y_rank1, n_pc = 1, center = FALSE)
  expect_lt(max(abs(resid)), 1e-10)
  expect_equal(attr(resid, "pc_removed"), 1L)
  expect_equal(rownames(resid), rownames(Y_rank1))
  expect_equal(colnames(resid), colnames(Y_rank1))
})

test_that("remove_EWVDS_pcs errors on non-finite input or n_pc < 1", {
  Y <- matrix(rnorm(20), nrow = 5, ncol = 4, dimnames = list(NULL, paste0("d", 1:4)))
  Y_bad <- Y; Y_bad[1, 1] <- NA
  expect_error(remove_EWVDS_pcs(Y_bad, n_pc = 1), "NA/NaN/Inf")
  expect_error(remove_EWVDS_pcs(Y, n_pc = 0), "n_pc must be >= 1")
})

# ---- fit_EWVDS_decomposition: core random-effects decomposition ----
# Mirrors the associated paper's Simulation 1 design (Methods/Simulation
# section): 20% signal edges with known true tau2, 80% null edges with
# tau2 = 0, so Storey's pi0 estimation (used unconditionally inside
# fit_EWVDS_decomposition for q_Q_storey_e) has enough null mass to
# converge -- an all-signal population makes qvalue::qvalue() itself
# error ("estimated pi0 <= 0"), independent of anything under test here.
make_mixed_fit_inputs <- function(seed = 42, E = 200, D = 20,
                                   tau2_true_signal = 0.03, frac_signal = 0.2) {
  set.seed(seed)
  donor_ncells <- setNames(rep(103, D), paste0("d", seq_len(D)))  # n_d - 3 = 100
  S2_mat <- matrix(1 / 100, nrow = E, ncol = D, dimnames = list(NULL, names(donor_ncells)))

  n_signal <- round(E * frac_signal)
  is_signal <- c(rep(FALSE, E - n_signal), rep(TRUE, n_signal))
  tau2_true_vec <- ifelse(is_signal, tau2_true_signal, 0)
  alpha_true <- rnorm(E, 0, 0.3)

  Y_mat <- matrix(NA_real_, E, D, dimnames = list(NULL, names(donor_ncells)))
  for (d in seq_len(D)) {
    Y_mat[, d] <- alpha_true + rnorm(E, 0, sqrt(tau2_true_vec)) + rnorm(E, 0, sqrt(1 / 100))
  }

  list(Y_mat = Y_mat, S2_mat = S2_mat, donor_ncells = donor_ncells,
       is_signal = is_signal, tau2_true_signal = tau2_true_signal)
}

test_that("fit_EWVDS_decomposition recovers near-zero tau2 on null edges and the known injected tau2 on signal edges", {
  inp <- make_mixed_fit_inputs()
  fit <- fit_EWVDS_decomposition(
    inp$Y_mat, inp$S2_mat, inp$donor_ncells,
    noise_method = "fisher_z", compute_Q0_sim = FALSE,
    evidence_p_method = "Q0_DL", excess_variance_method = "DL_tau",
    do_estimate_lfdr = FALSE
  )

  expect_equal(median(fit$tau2_DL_e[!inp$is_signal]), 0)
  expect_equal(median(fit$tau2_DL_e[inp$is_signal]), inp$tau2_true_signal, tolerance = 0.05)
})

test_that("fit_EWVDS_decomposition's floored quantities respect their stated bounds", {
  inp <- make_mixed_fit_inputs()
  fit <- fit_EWVDS_decomposition(
    inp$Y_mat, inp$S2_mat, inp$donor_ncells,
    noise_method = "fisher_z", compute_Q0_sim = FALSE,
    evidence_p_method = "Q0_DL", excess_variance_method = "DL_tau",
    do_estimate_lfdr = FALSE
  )

  expect_true(all(fit$phi2_e >= 0))
  expect_true(all(fit$tau2_e >= 0))
  expect_true(all(fit$Q0_e >= 0))
  expect_true(all(fit$lambda_DL_mat >= 0 & fit$lambda_DL_mat <= 1))
  expect_true(all(fit$p_Q_e >= 0 & fit$p_Q_e <= 1))
})

test_that("fit_EWVDS_decomposition errors rather than silently returning all-NA p_Q_e when evidence_p_method = 'Q0_sim' but compute_Q0_sim = FALSE", {
  inp <- make_mixed_fit_inputs()
  expect_error(
    fit_EWVDS_decomposition(
      inp$Y_mat, inp$S2_mat, inp$donor_ncells,
      compute_Q0_sim = FALSE, evidence_p_method = "Q0_sim",
      do_estimate_lfdr = FALSE
    ),
    "compute_Q0_sim = TRUE"
  )
})

test_that("fit_EWVDS_decomposition's closed-form DL path is deterministic given fixed inputs", {
  inp <- make_mixed_fit_inputs()
  args <- list(Y_mat = inp$Y_mat, S2_mat = inp$S2_mat, donor_ncells = inp$donor_ncells,
               noise_method = "fisher_z", compute_Q0_sim = FALSE,
               evidence_p_method = "Q0_DL", excess_variance_method = "DL_tau",
               do_estimate_lfdr = FALSE)
  fit1 <- do.call(fit_EWVDS_decomposition, args)
  fit2 <- do.call(fit_EWVDS_decomposition, args)
  expect_identical(fit1$tau2_DL_e, fit2$tau2_DL_e)
  expect_identical(fit1$Y_DL, fit2$Y_DL)
})

# ---- run_ewvds_stage / run_EWVDS_pipeline: small end-to-end integration ----

make_tiny_ewvds_fixture <- function(seed = 11, genes = c("g1", "g2", "g3", "g4"),
                                     donors = paste0("d", 1:6), n_cell_per_donor = 60) {
  set.seed(seed)
  comb_path <- tempfile("popcorn_test_")
  dir.create(comb_path)

  Sigma_true <- diag(length(genes)); Sigma_true[1, 2] <- Sigma_true[2, 1] <- 0.3
  L <- chol(Sigma_true)
  all_Z <- matrix(rnorm(length(donors) * n_cell_per_donor * length(genes)), ncol = length(genes)) %*% L
  colnames(all_Z) <- genes
  donor_id <- rep(donors, each = n_cell_per_donor)

  saveRDS(cor(all_Z), file.path(comb_path, "Sigma_global.rds"))
  write.csv(all_Z, file.path(comb_path, "Z_randomized.csv"), row.names = TRUE)

  dir.create(file.path(comb_path, "donors"))
  for (d in donors) {
    Zd <- all_Z[donor_id == d, , drop = FALSE]
    dd <- file.path(comb_path, "donors", d)
    dir.create(dd)
    write.csv(Zd, file.path(dd, "Z.csv"), row.names = TRUE)
    write.csv(cor(Zd), file.path(dd, "Sigma_Z.csv"), row.names = TRUE)
  }

  list(comb_path = comb_path, genes = genes, donors = donors)
}

test_that("run_ewvds_stage produces the expected output files and dimensions on a small fixture", {
  fx <- make_tiny_ewvds_fixture()
  on.exit(unlink(fx$comb_path, recursive = TRUE), add = TRUE)

  res <- run_ewvds_stage(comb_path = fx$comb_path, donors = fx$donors,
                          compute_Q0_sim = FALSE, do_estimate_lfdr = FALSE,
                          n_cores = 1, pipeline_tag = "unittest")

  expect_true(dir.exists(res$base_dir))
  expect_true(file.exists(file.path(res$base_dir, "EWVDS_result.rds")))
  for (d in fx$donors) {
    expect_true(file.exists(file.path(res$base_dir, "donor", d, "Sigma_DL.rds")))
  }

  result <- readRDS(file.path(res$base_dir, "EWVDS_result.rds"))
  n_edge <- choose(length(fx$genes), 2)
  expect_equal(dim(result$R_DL_mat), c(n_edge, length(fx$donors)))
  expect_true(all(result$R_DL_mat >= -1 & result$R_DL_mat <= 1))
})

test_that("run_ewvds_stage auto-detects donors from comb_path/donors/ when donors = NULL", {
  fx <- make_tiny_ewvds_fixture()
  on.exit(unlink(fx$comb_path, recursive = TRUE), add = TRUE)

  res <- run_ewvds_stage(comb_path = fx$comb_path, donors = NULL,
                          compute_Q0_sim = FALSE, do_estimate_lfdr = FALSE,
                          n_cores = 1, pipeline_tag = "unittest_autodetect")
  expect_equal(sort(res$inputs$donors), sort(fx$donors))
})

# ---- save_EWVDS_donor_networks: ensure_psd ----

# Three genes, two donors. Donor "bad" has edge values (0.9, 0.9, -0.9),
# which assemble into a non-PSD matrix; donor "good" is PSD.
make_psd_save_inputs <- function() {
  idx <- matrix(c(1, 1, 2,
                  2, 3, 3), ncol = 2)
  Y_DL <- cbind(bad = atanh(c(0.9, 0.9, -0.9)), good = atanh(c(0.1, 0.2, 0.1)))
  list(base_dir = tempfile("popcorn_psd_"), donors = c("bad", "good"), idx = idx,
       Z_g_edge = rep(0, 3), Y_DL = Y_DL, gene_names = paste0("g", 1:3))
}

call_save <- function(inp, ...) {
  save_EWVDS_donor_networks(base_dir = inp$base_dir, donors = inp$donors,
                            Sigma_global_Z = NULL, Z_g_edge = inp$Z_g_edge,
                            Y_shrunk = inp$Y_DL, idx = inp$idx,
                            gene_names = inp$gene_names, Y_DL = inp$Y_DL, ...)
}

test_that("save_EWVDS_donor_networks writes no PSD outputs by default", {
  inp <- make_psd_save_inputs()
  on.exit(unlink(inp$base_dir, recursive = TRUE), add = TRUE)
  call_save(inp)
  expect_true(file.exists(file.path(inp$base_dir, "donor", "bad", "Sigma_DL.rds")))
  expect_false(file.exists(file.path(inp$base_dir, "donor", "bad", "Sigma_DL_nearPD.rds")))
  expect_false(file.exists(file.path(inp$base_dir, "psd_check_DL.csv")))
})

test_that("save_EWVDS_donor_networks(ensure_psd = TRUE) projects only non-PSD donors and keeps Sigma_DL unchanged", {
  inp <- make_psd_save_inputs()
  on.exit(unlink(inp$base_dir, recursive = TRUE), add = TRUE)
  expect_message(call_save(inp, ensure_psd = TRUE), "1 of 2 donor matrices")

  dd <- function(d, f) readRDS(file.path(inp$base_dir, "donor", d, f))
  # Sigma_DL.rds still holds the raw edge-wise estimates
  expect_equal(dd("bad", "Sigma_DL.rds")["g2", "g3"], -0.9)
  # good donor: projected file is an exact copy
  expect_identical(dd("good", "Sigma_DL_nearPD.rds"), dd("good", "Sigma_DL.rds"))
  # bad donor: projected file is a valid correlation matrix
  P <- dd("bad", "Sigma_DL_nearPD.rds")
  expect_gte(min(eigen(P, only.values = TRUE)$values), -1e-8)
  expect_equal(unname(diag(P)), rep(1, 3))

  chk <- read.csv(file.path(inp$base_dir, "psd_check_DL.csv"))
  expect_equal(chk$donor, c("bad", "good"))
  expect_equal(chk$projected, c(TRUE, FALSE))
  expect_equal(chk$max_abs_change[2], 0)
})

test_that("run_ewvds_stage passes ensure_psd through to the saved donor networks", {
  fx <- make_tiny_ewvds_fixture()
  on.exit(unlink(fx$comb_path, recursive = TRUE), add = TRUE)
  res <- run_ewvds_stage(comb_path = fx$comb_path, donors = fx$donors,
                          compute_Q0_sim = FALSE, do_estimate_lfdr = FALSE,
                          n_cores = 1, pipeline_tag = "unittest_psd", ensure_psd = TRUE)
  expect_true(file.exists(file.path(res$base_dir, "psd_check_DL.csv")))
  for (d in fx$donors) {
    expect_true(file.exists(file.path(res$base_dir, "donor", d, "Sigma_DL_nearPD.rds")))
  }
})
