# ============================================================
# Stage 1b: Compute raw correlations (global and per-donor)
#
# Computes, for every edge, the global correlation from all cells and
# each donor's own correlation from that donor's cells alone (Methods,
# Figure 1b-IV). This is still part of Stage 1 -- it operates on the
# Gaussian-copula scores produced by run_marginal_stage() and does not
# yet involve the random-effects heterogeneity/shrinkage model of
# Stage 2.
# ============================================================

#' Compute global and per-donor raw correlations from Stage 1 scores
#'
#' Reads the Gaussian-copula scores written by [run_marginal_stage()],
#' computes the global correlation matrix from all cells pooled, then
#' splits cells by donor and computes each donor's own raw correlation
#' matrix from that donor's cells alone. Both use the same score matrix
#' `z_csv`, whose PIT values [run_marginal_stage()] has already clipped
#' to `[1e-6, 1 - 1e-6]`, so global and donor correlations share one
#' clipping threshold. This is Stage 1b: the bridge
#' between Stage 1's per-gene marginal fit and Stage 2's cross-donor
#' heterogeneity estimation, which reads the per-donor `Sigma_Z.csv`
#' files written here.
#'
#' @param comb_path Root output directory, same as passed to
#'   [run_marginal_stage()].
#' @param donor_col Column of the saved `metadata.rds$cell_covariates`
#'   holding donor identity. No default -- POPCORN estimates
#'   donor-*specific* correlations, so which column identifies donors
#'   is a modeling decision every caller must make explicitly.
#' @param u_csv Filename of the PIT uniform scores (relative to
#'   `comb_path`), as written by [run_marginal_stage()]. Only split into
#'   the per-donor `U.csv` files; not used for any correlation.
#' @param z_csv Filename of the Gaussian-copula scores (relative to
#'   `comb_path`), used for both the global and the per-donor
#'   correlations.
#' @param eps Deprecated and ignored. Earlier versions re-clipped `U` at
#'   this value before computing the global correlation only, so global
#'   and donor correlations used different clipping thresholds. Passing
#'   any value now gives a warning.
#' @param n_cores Number of parallel workers for the per-donor split.
#'
#' @return Invisibly, a list with `comb_path`, `donors`, and
#'   `donor_outdir`.
#'
#' @section Output:
#' Writes `comb_path/Sigma_global.rds` (genes x genes, skipped if it
#' already exists) and, for every donor, `comb_path/donors/<donor_id>/`
#' containing `U.csv`, `Z.csv`, and `Sigma_Z.csv` (that donor's own
#' genes x genes raw correlation matrix).
#'
#' @export
run_donor_split <- function(
  comb_path,
  donor_col,
  u_csv     = "U_randomized.csv",
  z_csv     = "Z_randomized.csv",
  eps       = NULL,
  n_cores   = 40L
) {
  if (!is.null(eps)) {
    warning("run_donor_split(): `eps` is deprecated and ignored. Global and ",
            "donor correlations are both computed from `z_csv`, whose PIT ",
            "values run_marginal_stage() clips to [1e-6, 1 - 1e-6].",
            call. = FALSE)
  }
  U        <- utils::read.csv(file.path(comb_path, u_csv), row.names = 1, check.names = FALSE)
  Z        <- utils::read.csv(file.path(comb_path, z_csv), row.names = 1, check.names = FALSE)
  meta     <- readRDS(file.path(comb_path, "metadata.rds"))
  cell_cov <- meta$cell_covariates

  # --- global Sigma ---
  # Same Z matrix as the per-donor correlations below, so global and donor
  # correlations share one clipping threshold (Methods, Stage 1).
  sig_path <- file.path(comb_path, "Sigma_global.rds")
  if (!file.exists(sig_path)) {
    Sigma_global <- cor(as.matrix(Z), use = "pairwise.complete.obs")
    diag(Sigma_global) <- 0
    saveRDS(Sigma_global, sig_path)
    message(sprintf("[Stage 1: raw correlation] Saved Sigma_global.rds  (%d NAs in off-diagonal)",
                    sum(is.na(Sigma_global))))
  } else {
    message("[Stage 1: raw correlation] Sigma_global.rds already exists, skipping.")
  }

  # --- per-donor split ---
  donor_outdir <- file.path(comb_path, "donors")
  dir.create(donor_outdir, showWarnings = FALSE, recursive = TRUE)
  donors <- unique(cell_cov[[donor_col]])
  message(sprintf("[Stage 1: raw correlation] Splitting %d donors ...", length(donors)))

  process_one_donor <- function(donor_id) {
    idx <- which(cell_cov[[donor_col]] == donor_id)
    if (length(idx) == 0L) {
      warning("No cells for donor: ", donor_id)
      return(NULL)
    }
    U_d     <- U[idx, , drop = FALSE]
    Z_d     <- Z[idx, , drop = FALSE]
    Sigma_d <- cor(Z_d, use = "pairwise.complete.obs")

    donor_dir <- file.path(donor_outdir, donor_id)
    dir.create(donor_dir, showWarnings = FALSE, recursive = TRUE)
    utils::write.csv(U_d,     file.path(donor_dir, "U.csv"),       row.names = TRUE, quote = FALSE)
    utils::write.csv(Z_d,     file.path(donor_dir, "Z.csv"),       row.names = TRUE, quote = FALSE)
    utils::write.csv(Sigma_d, file.path(donor_dir, "Sigma_Z.csv"), row.names = TRUE, quote = FALSE)
    invisible(TRUE)
  }

  parallel::mclapply(donors, process_one_donor, mc.cores = n_cores)

  message("[Stage 1: raw correlation] Done.")
  invisible(list(comb_path = comb_path, donors = donors, donor_outdir = donor_outdir))
}
