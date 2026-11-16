# ============================================================
# Complete pipeline: Stage 1 (marginal fit + raw correlation)
# + Stage 2 (donor heterogeneity estimation and shrinkage)
# ============================================================

#' Run the complete POPCORN pipeline
#'
#' Top-level orchestrator matching the two-stage structure described in
#' the paper's Methods section: Stage 1 (marginal fitting and raw
#' correlation, via [run_marginal_stage()] and [run_donor_split()]) then
#' Stage 2 (donor heterogeneity estimation and shrinkage, via
#' [run_ewvds_stage()]).
#'
#' @param sce A `SingleCellExperiment`, genes x cells, with all of
#'   `covariates` present in `colData(sce)`. Should already be subset
#'   to a single cell type: POPCORN's donor-level correlation
#'   structure is defined within one cell type, and `celltype` below
#'   does not perform that subsetting for you.
#' @param comb_path Root output directory for every stage's intermediate
#'   and final output (created if absent).
#' @param donor_col `colData(sce)` column holding donor identity. No
#'   default -- POPCORN estimates donor-*specific* correlations, so
#'   which column identifies donors is a modeling decision every
#'   caller must make explicitly. We recommend also including it in
#'   `covariates`/`mu_formula`; see [run_marginal_stage()].
#' @param gene_names Character vector of genes to analyze; must be a
#'   subset of `rownames(sce)`. If `NULL` (default), every gene in
#'   `rownames(sce)` is used.
#' @param covariates Character vector of `colData(sce)` columns to pass
#'   to the Stage 1 marginal model. Default `"log1p_total_counts"`
#'   alone; we recommend also adding `donor_col`; see
#'   [run_marginal_stage()].
#' @param mu_formula Stage 1 GAMLSS mu formula (right-hand side only).
#'   Default `"log1p_total_counts"` alone, matching the default
#'   `covariates`.
#' @param sigma_formula Stage 1 GAMLSS sigma formula. Default `"1"`.
#' @param family_use Stage 1 marginal family: `"nb"`, `"poisson"`,
#'   `"gaussian"`, or `"auto"`.
#' @param celltype `colData(sce)` column holding a cell-type-like
#'   label; see [run_marginal_stage()] for details and the default.
#' @param n_cores_marginal Parallel workers for Stage 1's per-gene fit.
#' @param timeout_sec Per-gene timeout, in seconds, for Stage 1.
#' @param blas_threads_per_gene Caps BLAS/OpenMP threads per Stage 1 gene
#'   fit (see [run_marginal_stage()]).
#' @param u_csv,z_csv Filenames (relative to `comb_path`) of Stage 1's
#'   PIT uniform / Gaussian-copula score matrices.
#' @param eps_pit Deprecated and ignored (see [run_donor_split()]). Global
#'   and donor correlations both use Stage 1's scores, whose PIT values
#'   are clipped to `[1e-6, 1 - 1e-6]`. Passing any value gives a warning.
#' @param n_cores_split Parallel workers for Stage 1b's per-donor split.
#' @param donors Character vector of donor IDs for Stage 2; if `NULL`,
#'   auto-detected from `comb_path/donors/`.
#' @param noise_method Stage 2 per-donor sampling-variance estimator:
#'   `"fisher_z"` (closed-form, used throughout the associated paper) or
#'   `"bootstrap"` (alternative; not used to produce any result in that
#'   paper and not independently tested -- see
#'   `bootstrap_edge_variance()`).
#' @param bootstrap_weight_method Passed through to the bootstrap
#'   variance path when `noise_method = "bootstrap"`.
#' @param do_pc_removal,n_pc If `do_pc_removal = TRUE`, removes the
#'   leading `n_pc` principal components from the donor-deviation matrix
#'   before the heterogeneity decomposition (see
#'   `remove_EWVDS_pcs()`). This is an exploratory option: it is not
#'   described in the associated paper and no analysis there sets it to
#'   `TRUE`.
#' @param seed_bootstrap Random seed for the bootstrap variance path.
#' @param compute_Q0_sim Whether to additionally compute a
#'   simulation-based (pseudo-donor) null distribution for `Q_e`.
#' @param null_method `"pseudo_donor"` (default) or `"gaussian"`,
#'   controlling how `compute_Q0_sim` constructs its null.
#' @param pseudo_donor_replace,pseudo_donor_mode Control how pseudo-donor
#'   bucket membership is drawn. `pseudo_donor_mode = "permutation"`
#'   (default, used in the associated paper) partitions cells into
#'   disjoint pseudo-donors. `"resampling"` draws each pseudo-donor's
#'   cells independently (with or without replacement, per
#'   `pseudo_donor_replace`), which biases `mean(Q_e^null)` slightly
#'   below `D - 1`. Only used when `compute_Q0_sim = TRUE`.
#' @param n_null_sim Number of pseudo-donor null replicates.
#' @param seed_Q0_sim Random seed for the pseudo-donor null simulation.
#' @param evidence_p_method Which `Q_e` null to use for the reported
#'   edge-level p-value: `"Q0_DL"` (asymptotic chi-squared) or `"Q0_sim"`
#'   (the simulation-based null above).
#' @param excess_variance_method Which excess-variance estimator feeds
#'   the shrinkage weight: `"DL_tau"`, `"DL_phi"`, `"sim_tau"`, or
#'   `"sim_phi"`.
#' @param phi_alpha Shrinkage-weight tuning parameter.
#' @param do_estimate_lfdr Whether to additionally estimate a
#'   kernel-smoothed local false discovery rate.
#' @param lfdr_pi0_method,lfdr_pi0_lambda,lfdr_K,lfdr_kernel Local-fdr
#'   estimation settings (see `estimate_lfdr()`).
#' @param min_var,min_C,eps_ewvds Numerical floors used in Stage 2's
#'   random-effects decomposition.
#' @param pipeline_tag String appended to the Stage 2 output subdirectory
#'   name, to distinguish parameter configurations run against the same
#'   `comb_path`.
#' @param n_cores_ewvds Parallel workers for Stage 2.
#' @param force_recompute_stage2 If `FALSE` (default), reuse cached
#'   `Sigma_donor.rds`/`bootstrap_var_fisherZ.rds`/`Z_randomized.rds`
#'   under `comb_path/donor_network_EWVDS*/` where present, instead of
#'   recomputing them. This is Stage 2 only -- Stage 1a's per-gene fits
#'   and Stage 1b's `Sigma_global.rds`/per-donor files are always cached
#'   by filename with no override switch; `run_marginal_stage()` errors
#'   instead if `comb_path` already has Stage 1 output for a different
#'   `sce` (see its own docs).
#' @param ensure_psd If `TRUE`, also write a positive semi-definite
#'   version of each donor's correlation matrix, `Sigma_DL_nearPD.rds`
#'   (see [save_EWVDS_donor_networks()] and [ensure_psd_correlation()]).
#'   Edge-level estimates are unaffected. Default `FALSE`.
#' @param fix_rng_streams Passed to [run_marginal_stage()] (see that
#'   function's documentation). Keep the default `TRUE`.
#' @param seed Optional integer, passed to [run_marginal_stage()]. Seeds
#'   Stage 1's randomized PIT so the whole pipeline is bit-reproducible,
#'   including across repeated calls in one R session. Stage 2's
#'   pseudo-donor null (`compute_Q0_sim = TRUE`) is seeded separately by
#'   `seed_Q0_sim`. Default `NULL` leaves the RNG untouched.
#'
#' @return Invisibly, Stage 2's result object (see
#'   `run_EWVDS_pipeline()`).
#' @export
run_coexp_pipeline <- function(
  # --- Stage 1a: marginal fit ---
  sce,
  comb_path,
  donor_col,
  gene_names    = NULL,
  covariates    = c("log1p_total_counts"),
  mu_formula    = "log1p_total_counts",
  sigma_formula = "1",
  family_use       = "nb",
  celltype         = "cell_type",
  n_cores_marginal = 20L,
  timeout_sec      = 3600L,
  blas_threads_per_gene = 2L,

  # --- Stage 1b: raw correlation (donor split) ---
  u_csv            = "U_randomized.csv",
  z_csv            = "Z_randomized.csv",
  eps_pit          = NULL,
  n_cores_split    = 40L,

  # --- Stage 2: heterogeneity estimation and shrinkage ---
  donors                  = NULL,
  noise_method            = "fisher_z",
  bootstrap_weight_method = "fisher_z",
  do_pc_removal           = FALSE,
  n_pc                    = 5,
  seed_bootstrap          = 42,
  compute_Q0_sim          = FALSE,
  null_method             = "pseudo_donor",
  pseudo_donor_replace    = FALSE,
  pseudo_donor_mode       = "permutation",
  n_null_sim              = 999,
  seed_Q0_sim             = 43,
  evidence_p_method       = "Q0_DL",
  excess_variance_method  = "DL_tau",
  phi_alpha               = 0.5,
  do_estimate_lfdr        = TRUE,
  lfdr_pi0_method         = "storey",
  lfdr_pi0_lambda         = 0.5,
  lfdr_K                  = 31,
  lfdr_kernel             = "triangular",
  min_var                 = 1e-10,
  min_C                   = 1e-12,
  eps_ewvds               = 1e-12,
  pipeline_tag            = "pipeline_DL_noisefisher",
  n_cores_ewvds           = 30L,
  force_recompute_stage2         = FALSE,
  ensure_psd              = FALSE,
  fix_rng_streams         = TRUE,  # see run_marginal_stage()
  seed                    = NULL   # see run_marginal_stage(): exact reproducibility of Stage 1's PIT
) {
  if (!is.null(eps_pit)) {
    warning("run_coexp_pipeline(): `eps_pit` is deprecated and ignored. Global ",
            "and donor correlations both use Stage 1's scores, whose PIT values ",
            "are clipped to [1e-6, 1 - 1e-6].", call. = FALSE)
  }
  message("========== run_coexp_pipeline START ==========")
  message("comb_path: ", comb_path)

  message("\n--- Stage 1a: Marginal fit (scDesign3) ---")
  run_marginal_stage(
    sce           = sce,
    comb_path     = comb_path,
    gene_names    = gene_names,
    covariates    = covariates,
    mu_formula    = mu_formula,
    sigma_formula = sigma_formula,
    family_use    = family_use,
    celltype      = celltype,
    donor_col     = donor_col,
    n_cores       = n_cores_marginal,
    timeout_sec   = timeout_sec,
    blas_threads_per_gene = blas_threads_per_gene,
    fix_rng_streams = fix_rng_streams,
    seed          = seed
  )

  message("\n--- Stage 1b: Raw correlation (global + per-donor) ---")
  run_donor_split(
    comb_path = comb_path,
    donor_col = donor_col,
    u_csv     = u_csv,
    z_csv     = z_csv,
    n_cores   = n_cores_split
  )
  saveRDS(Sys.time(), file.path(comb_path, "time_stage1b_done.rds"))

  message("\n--- Stage 2: Donor heterogeneity estimation and shrinkage ---")
  res <- run_ewvds_stage(
    comb_path               = comb_path,
    donors                  = donors,
    noise_method            = noise_method,
    bootstrap_weight_method = bootstrap_weight_method,
    do_pc_removal           = do_pc_removal,
    n_pc                    = n_pc,
    seed_bootstrap          = seed_bootstrap,
    compute_Q0_sim          = compute_Q0_sim,
    null_method             = null_method,
    pseudo_donor_replace    = pseudo_donor_replace,
    pseudo_donor_mode       = pseudo_donor_mode,
    n_null_sim              = n_null_sim,
    seed_Q0_sim             = seed_Q0_sim,
    evidence_p_method       = evidence_p_method,
    excess_variance_method  = excess_variance_method,
    phi_alpha               = phi_alpha,
    do_estimate_lfdr        = do_estimate_lfdr,
    lfdr_pi0_method         = lfdr_pi0_method,
    lfdr_pi0_lambda         = lfdr_pi0_lambda,
    lfdr_K                  = lfdr_K,
    lfdr_kernel             = lfdr_kernel,
    min_var                 = min_var,
    min_C                   = min_C,
    eps                     = eps_ewvds,
    pipeline_tag            = pipeline_tag,
    n_cores                 = n_cores_ewvds,
    force_recompute_stage2         = force_recompute_stage2,
    ensure_psd              = ensure_psd
  )

  message("========== run_coexp_pipeline DONE ==========")
  invisible(res)
}
