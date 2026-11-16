# -----------------------------
# Loading Stage 1's per-donor outputs into Stage 2's edge x donor
# matrices, and saving Stage 2's shrunk per-donor correlation matrices.
# -----------------------------

#' Load Stage 1 outputs into Stage 2's edge x donor matrices
#'
#' Reads the global correlation matrix and every donor's raw correlation
#' matrix (written by [run_donor_split()]), Fisher-z transforms them,
#' and assembles the edge x donor deviation matrix `Y_mat` (each donor's
#' deviation from the global correlation, in Fisher-z space) and its
#' sampling-variance matrix `S2_mat` that Stage 2's heterogeneity
#' decomposition operates on.
#'
#' @param COMB_PATH Root output directory, same as passed to
#'   [run_marginal_stage()] / [run_donor_split()].
#' @param donors Character vector of donor IDs to load.
#' @param n_boot Number of bootstrap replicates, only used when
#'   `noise_method = "bootstrap"`.
#' @param noise_method `"fisher_z"` (closed-form sampling variance,
#'   `1 / (n_d - 3)`) or `"bootstrap"` (see `bootstrap_edge_variance()`;
#'   not used in the associated paper).
#' @param force_recompute_stage2 If `TRUE`, ignore any cached
#'   `Sigma_donor.rds`/`bootstrap_var_fisherZ.rds`/`Z_randomized.rds` and
#'   recompute from the underlying CSVs.
#' @param use_pairwise If `TRUE`, correlations use
#'   `"pairwise.complete.obs"` instead of `"everything"`.
#' @param eps Clip value used before the Fisher-z transform.
#' @param seed_bootstrap Base random seed for the bootstrap variance
#'   path (offset per donor).
#' @param n_cores Number of parallel workers across donors.
#' @param output_subdir Subdirectory of `COMB_PATH` to read/write
#'   Stage 2 outputs under.
#'
#' @return A list with (among others) `base_dir`, `Sigma_global`,
#'   `Sigma_global_Z`, `Z_g_edge`, `idx`, `gene_names`, `donors`,
#'   `donor_ncells`, `Y_mat`, `S2_mat`, `Z_cell_mat`.
#' @export
load_EWVDS_inputs <- function(
  COMB_PATH,
  donors,
  n_boot = 50,
  noise_method = c("fisher_z", "bootstrap"),
  force_recompute_stage2 = FALSE,
  use_pairwise = FALSE,
  eps = 1e-6,
  seed_bootstrap = 42,
  n_cores = 1,
  output_subdir = "donor_network_EWVDS"
) {
  noise_method <- match.arg(noise_method)

  base_dir <- file.path(COMB_PATH, output_subdir)
  dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)

  Sigma_global <- readRDS(file.path(COMB_PATH, "Sigma_global.rds"))
  gene_names <- colnames(Sigma_global)

  idx <- get_upper_idx(Sigma_global)
  E <- nrow(idx)
  D <- length(donors)

  Sigma_global_Z <- fisher_z(Sigma_global, eps = eps)
  Z_g_edge <- Sigma_global_Z[cbind(idx[, 1], idx[, 2])]

  # Process donors in parallel (donor-level parallelism avoids mclapply overhead
  # of per-bootstrap parallelism; bootstrap stays sequential within each worker).
  donor_results <- parallel::mclapply(seq_along(donors), function(i) {
    d <- donors[i]

    donor_raw_dir <- file.path(COMB_PATH, "donors", d)
    donor_out_dir <- file.path(base_dir, "donor", d)
    dir.create(donor_out_dir, recursive = TRUE, showWarnings = FALSE)

    sigma_file    <- file.path(donor_out_dir, "Sigma_donor.rds")
    boot_var_file <- file.path(donor_out_dir, "bootstrap_var_fisherZ.rds")

    Sigma_d_raw <- as.matrix(
      read.csv(file.path(donor_raw_dir, "Sigma_Z.csv"), row.names = 1, check.names = FALSE)
    )
    Z_donor <- as.matrix(
      read.csv(file.path(donor_raw_dir, "Z.csv"), row.names = 1, check.names = FALSE)
    )

    if (!setequal(rownames(Sigma_d_raw), gene_names) || !setequal(colnames(Sigma_d_raw), gene_names))
      stop(paste0("Donor ", d, ": Sigma_Z gene names do not match Sigma_global gene_names."))
    Sigma_d_raw <- Sigma_d_raw[gene_names, gene_names, drop = FALSE]

    if (!setequal(colnames(Z_donor), gene_names))
      stop(paste0("Donor ", d, ": Z.csv gene names do not match Sigma_global gene_names."))
    Z_donor <- Z_donor[, gene_names, drop = FALSE]

    donor_ncell <- nrow(Z_donor)

    if (!file.exists(sigma_file) || force_recompute_stage2) {
      saveRDS(Sigma_d_raw, sigma_file)
      Sigma_d <- Sigma_d_raw
    } else {
      Sigma_d <- readRDS(sigma_file)
    }

    Sigma_d_Z <- fisher_z(Sigma_d, eps = eps)
    Y_col <- Sigma_d_Z[cbind(idx[, 1], idx[, 2])] - Z_g_edge

    if (noise_method == "fisher_z") {
      if (donor_ncell <= 3)
        stop(paste0("Donor ", d, ": must have more than 3 cells for Fisher-Z variance."))
      s2_col <- rep(1 / (donor_ncell - 3), E)
    } else {
      if (!file.exists(boot_var_file) || force_recompute_stage2) {
        s2_col <- bootstrap_edge_variance(
          Z_donor      = Z_donor,
          idx          = idx,
          n_boot       = n_boot,
          seed         = seed_bootstrap + i,
          use_pairwise = use_pairwise,
          eps          = eps
        )
        saveRDS(s2_col, boot_var_file)
      } else {
        s2_col <- readRDS(boot_var_file)
      }
    }

    list(d = d, Y_col = Y_col, s2_col = s2_col,
         donor_ncell = donor_ncell, Sigma_d = Sigma_d)
  }, mc.cores = n_cores)

  failed <- vapply(donor_results, inherits, logical(1), what = "try-error")
  if (any(failed)) {
    first_err <- conditionMessage(attr(donor_results[[which(failed)[1]]], "condition"))
    stop(sprintf("load_EWVDS_inputs: %d / %d donor(s) failed.\nFirst error: %s",
                 sum(failed), D, first_err))
  }

  Y_mat          <- matrix(NA_real_, nrow = E, ncol = D, dimnames = list(NULL, donors))
  S2_mat         <- matrix(NA_real_, nrow = E, ncol = D, dimnames = list(NULL, donors))
  donor_ncells   <- setNames(rep(NA_integer_, D), donors)
  donor_sigma_list <- list()
  donor_var_list   <- list()

  for (i in seq_along(donors)) {
    d <- donors[i]
    Y_mat[, i]           <- donor_results[[i]]$Y_col
    S2_mat[, i]          <- donor_results[[i]]$s2_col
    donor_ncells[d]      <- donor_results[[i]]$donor_ncell
    donor_sigma_list[[d]]<- donor_results[[i]]$Sigma_d
    donor_var_list[[d]]  <- donor_results[[i]]$s2_col
  }

  # Read Z_cell_mat after the donor loop so it is not copied into fork memory.
  # Cache as RDS to avoid repeated slow CSV reads from network storage.
  z_cache_file <- file.path(base_dir, "Z_randomized.rds")
  if (!file.exists(z_cache_file) || force_recompute_stage2) {
    Z_cell_mat <- as.matrix(
      read.csv(file.path(COMB_PATH, "Z_randomized.csv"),
               row.names = 1, check.names = FALSE)
    )
    saveRDS(Z_cell_mat, z_cache_file)
  } else {
    Z_cell_mat <- readRDS(z_cache_file)
  }

  if (!identical(colnames(Z_cell_mat), gene_names))
    stop("colnames(Z_cell_mat) do not match gene_names from Sigma_global.")

  saveRDS(
    list(
      Sigma_global = Sigma_global, Sigma_global_Z = Sigma_global_Z,
      Z_g_edge = Z_g_edge, idx = idx, gene_names = gene_names,
      donors = donors, donor_ncells = donor_ncells,
      Y_mat = Y_mat, S2_mat = S2_mat,
      Z_cell_mat = Z_cell_mat, noise_method = noise_method
    ),
    file.path(base_dir, "EWVDS_inputs.rds")
  )

  list(
    base_dir = base_dir,
    Sigma_global = Sigma_global, Sigma_global_Z = Sigma_global_Z,
    Z_g_edge = Z_g_edge, idx = idx, gene_names = gene_names,
    donors = donors, donor_ncells = donor_ncells,
    Y_mat = Y_mat, S2_mat = S2_mat,
    Z_cell_mat = Z_cell_mat,
    donor_sigma_list = donor_sigma_list, donor_var_list = donor_var_list,
    noise_method = noise_method
  )
}

#' Save Stage 2's shrunk per-donor correlation matrices
#'
#' For every donor, reconstructs the full gene x gene shrunk correlation
#' matrix from the edge-level deviation `Y_shrunk` (and, if supplied,
#' `Y_DL`), adds back the global correlation, and inverse-Fisher-z
#' transforms to the correlation scale.
#'
#' @param base_dir Stage 2 output directory (`inputs$base_dir` from
#'   [load_EWVDS_inputs()]).
#' @param donors Character vector of donor IDs.
#' @param Sigma_global_Z Global correlation matrix, Fisher-z scale.
#' @param Z_g_edge Edge-level global correlation, Fisher-z scale
#'   (`Sigma_global_Z` indexed by `idx`).
#' @param Y_shrunk Edge x donor deviation matrix, the local-fdr-weighted
#'   (`rho_e`) shrinkage variant -- see
#'   `fit_EWVDS_decomposition()`'s output field glossary.
#' @param idx Two-column edge index matrix (see [get_upper_idx()]).
#' @param gene_names Optional character vector of gene names for the
#'   output matrices' row/column names.
#' @param Y_DL Optional edge x donor deviation matrix, the canonical
#'   POPCORN shrinkage variant (no `rho_e` weighting). If supplied, also
#'   writes `Sigma_DL.rds` per donor.
#' @param ensure_psd If `TRUE` (and `Y_DL` is supplied), also write a
#'   positive semi-definite version of each donor's `Sigma_DL` (see
#'   [ensure_psd_correlation()]). Edge-wise shrinkage does not guarantee
#'   that the assembled matrix is PSD, which matters only for downstream
#'   uses that need a valid correlation matrix (e.g. inversion or
#'   simulation). Default `FALSE`.
#' @param psd_tol Tolerance passed to [ensure_psd_correlation()].
#'
#' @return Invisible `TRUE`.
#'
#' @section Output:
#' Writes, per donor, `base_dir/donor/<donor_id>/Sigma_EWVDS.rds` and
#' (when `Y_DL` is supplied) `Sigma_DL.rds` -- the latter is the
#' canonical POPCORN per-donor correlation matrix. `Sigma_DL.rds` always
#' holds the edge-wise estimates exactly as in `EWVDS_result.rds`.
#'
#' With `ensure_psd = TRUE`, additionally writes per donor
#' `Sigma_DL_nearPD.rds` (identical to `Sigma_DL.rds` when that is
#' already PSD, the nearest correlation matrix otherwise), plus a
#' per-donor summary `base_dir/psd_check_DL.csv` (smallest eigenvalue
#' before/after, whether projected, largest absolute change). The legacy
#' `Sigma_EWVDS.rds` variant is never projected.
#' @export
save_EWVDS_donor_networks <- function(
  base_dir,
  donors,
  Sigma_global_Z,
  Z_g_edge,
  Y_shrunk,
  idx,
  gene_names = NULL,
  Y_DL = NULL,
  ensure_psd = FALSE,
  psd_tol = 1e-8
) {
  D <- length(donors)
  do_psd <- ensure_psd && !is.null(Y_DL)
  if (ensure_psd && is.null(Y_DL)) {
    warning("ensure_psd = TRUE has no effect without Y_DL.")
  }
  psd_rows <- vector("list", D)

  for (i in seq_len(D)) {
    d <- donors[i]
    donor_out_dir <- file.path(base_dir, "donor", d)
    dir.create(donor_out_dir, recursive = TRUE, showWarnings = FALSE)

    # EWS version (default)
    Z_EWVDS_edge <- Z_g_edge + Y_shrunk[, i]
    Sigma_EWVDS_Z <- edge_to_symmetric_mat(
      edge_values = Z_EWVDS_edge,
      idx = idx,
      gene_names = gene_names,
      diag_value = 0
    )
    Sigma_EWVDS <- inv_fisher_z(Sigma_EWVDS_Z)
    diag(Sigma_EWVDS) <- 1
    saveRDS(Sigma_EWVDS, file.path(donor_out_dir, "Sigma_EWVDS.rds"))

    # DL version (optional)
    if (!is.null(Y_DL)) {
      Z_DL_edge <- Z_g_edge + Y_DL[, i]
      Sigma_DL_Z <- edge_to_symmetric_mat(
        edge_values = Z_DL_edge,
        idx = idx,
        gene_names = gene_names,
        diag_value = 0
      )
      Sigma_DL <- inv_fisher_z(Sigma_DL_Z)
      diag(Sigma_DL) <- 1
      saveRDS(Sigma_DL, file.path(donor_out_dir, "Sigma_DL.rds"))

      if (do_psd) {
        chk <- ensure_psd_correlation(Sigma_DL, tol = psd_tol)
        saveRDS(chk$mat, file.path(donor_out_dir, "Sigma_DL_nearPD.rds"))
        psd_rows[[i]] <- data.frame(
          donor = d, projected = chk$projected,
          min_eig_before = chk$min_eig_before, min_eig_after = chk$min_eig_after,
          max_abs_change = chk$max_abs_change
        )
      }
    }
  }

  if (do_psd) {
    psd_df <- do.call(rbind, psd_rows)
    utils::write.csv(psd_df, file.path(base_dir, "psd_check_DL.csv"), row.names = FALSE)
    if (any(psd_df$projected)) {
      message(sprintf(
        "[Stage 2] %d of %d donor matrices were not PSD and were projected (largest entry change %.2g); see psd_check_DL.csv.",
        sum(psd_df$projected), D, max(psd_df$max_abs_change)))
    }
  }

  invisible(TRUE)
}
