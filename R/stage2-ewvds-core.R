# ============================================================
# EWVDS Pipeline
# Edge-wise Variance Decomposition Shrinkage
#
# Stage 2's core random-effects decomposition: separates genuine
# between-donor heterogeneity from sampling noise at each edge (the
# DerSimonian-Laird-style Q_e statistic and tau2_e/phi2_e heterogeneity
# estimate), then adaptively shrinks each donor's raw deviation from the
# global correlation toward 0 by a ratio set by that heterogeneity
# relative to the donor's own sampling variance.
# ============================================================

# ============================================================
# OUTPUT FIELD GLOSSARY (read this before using any field returned by
# fit_EWVDS_decomposition() or saved in EWVDS_result.rds /
# EWVDS_full_intermediate.rds -- several of these fields are easy to
# misread as something they are not):
#
# Edge-local, correlation-scale (safe to compare a gene/edge subset's
# output against a full-panel run):
#   Sigma_global / Sigma_global_Z -- the global (all-cells-pooled,
#     pre-shrinkage) raw correlation / its Fisher-z transform. THIS is
#     "the raw correlation," if you want one number per edge.
#   Z_g_edge -- Sigma_global_Z sliced to one edge (same thing, indexed).
#
# NOT edge-local, do not compare at reduced gene-panel size (these pool
# information across ALL edges via density estimation / the full null
# distribution):
#   lfdr_e, q_Q_BH_e, q_Q_storey_e, p_Q_e, lfdr_fit, pi0_Q_storey.
#
# DEVIATIONS, not correlations (per-donor minus Z_g_edge, in Fisher-z
# space -- small numbers clustered near 0 even for real signal edges; do
# not read these as "the per-donor correlation"):
#   Y_mat    -- raw per-donor deviation (one column per donor).
#   Y_shrunk -- Y_mat after reliability-weighted shrinkage, ADDITIONALLY
#               scaled by rho_e (= 1 - lfdr_e, see below) -- an
#               alternative variant, not the one this package calls
#               "POPCORN."
#   Y_DL     -- Y_mat after reliability-weighted shrinkage ONLY, no
#               rho_e/lfdr_e weighting -- THIS is the canonical POPCORN
#               per-donor estimate (pre-correlation-scale).
#
# Absolute, correlation-scale, SHRUNK per-donor estimates (add Z_g_edge
# back onto the deviation, then inverse-Fisher-z; see
# Z_EWVDS_mat/R_EWVDS_mat/Z_DL_mat/R_DL_mat below) -- use R_DL_mat
# unless you specifically want the lfdr-weighted variant:
#   Z_DL_mat -- Y_DL + Z_g_edge (Fisher-z scale).
#   R_DL_mat -- inv_fisher_z(Z_DL_mat) -- the canonical POPCORN
#               per-donor correlation estimate, correlation scale.
#   Z_EWVDS_mat / R_EWVDS_mat -- same, from Y_shrunk instead of Y_DL.
#
# rho_e -- NOT a correlation of any kind. It is 1 - p_evidence_e (an
#   lfdr-derived evidence weight), used only to build the
#   Y_shrunk/R_EWVDS_mat variant above. It is NOT edge-local (depends on
#   the lfdr fit across all edges) and should never be read as "the
#   raw" or "the shrunk" correlation.
# ============================================================

#' Fit Stage 2's heterogeneity decomposition and shrinkage for every
#' edge
#'
#' The core random-effects computation: a DerSimonian-Laird-style
#' weighted-mean/Q-statistic step (`alpha0_e`, `Q0_e`, `tau2_DL_e`),
#' optionally refined against a simulation-based (pseudo-donor) null
#' instead of the asymptotic chi-squared null, followed by a one-step
#' re-weighted geometry (`alpha1_e`, `Q1_e`) and the final
#' reliability-weighted shrinkage producing `Y_DL` (the canonical
#' POPCORN deviation estimate) and `Y_shrunk` (the local-fdr-weighted
#' variant). See the file-level output field glossary above before
#' reading any returned field.
#'
#' @param Y_mat Edge x donor deviation matrix, Fisher-z scale (from
#'   [load_EWVDS_inputs()]).
#' @param S2_mat Edge x donor sampling-variance matrix, same shape as
#'   `Y_mat`.
#' @param donor_ncells Named integer vector of per-donor cell counts,
#'   matching `colnames(Y_mat)`.
#' @param noise_method `"fisher_z"` (closed-form, used throughout the
#'   associated paper) or `"bootstrap"` (see [bootstrap_edge_variance()];
#'   not used in that paper and not independently tested).
#' @param bootstrap_weight_method Only used when `noise_method` or the
#'   weighting scheme requests a bootstrap-derived weight.
#' @param compute_Q0_sim Whether to additionally compute a
#'   simulation-based (pseudo-donor) null for `Q_e`. Must be `TRUE` if
#'   `evidence_p_method = "Q0_sim"` -- errors otherwise, rather than
#'   silently leaving `p_Q_e` all `NA`.
#' @param null_method `"pseudo_donor"` (default) or `"gaussian"`,
#'   controlling how `compute_Q0_sim` constructs its null.
#' @param pseudo_donor_mode,pseudo_donor_replace Passed to
#'   [generate_pseudodonor_indices()] when `null_method = "pseudo_donor"`.
#'   Default mode `"permutation"` (disjoint pseudo-donors, used in the
#'   associated paper); see [generate_pseudodonor_indices()].
#' @param Z_cell_mat,Z_g_edge,idx Required when
#'   `null_method = "pseudo_donor"`; see [load_EWVDS_inputs()].
#' @param use_pairwise If `TRUE`, correlations in the pseudo-donor null
#'   use `"pairwise.complete.obs"` instead of `"everything"`.
#' @param n_cores Parallel workers for the pseudo-donor null simulation.
#' @param n_null_sim Number of pseudo-donor null replicates.
#' @param null_seed Random seed for the null simulation.
#' @param evidence_p_method Which `Q_e` null to use for the reported
#'   edge-level p-value: `"Q0_DL"` (asymptotic chi-squared, used
#'   throughout the associated paper) or `"Q0_sim"` (the simulation-based
#'   null above).
#' @param excess_variance_method Which excess-variance estimator feeds
#'   the shrinkage weight: `"DL_tau"` (used throughout the associated
#'   paper), `"DL_phi"`, `"sim_tau"`, or `"sim_phi"`.
#' @param phi_alpha Shrinkage-weight tuning parameter (quantile level
#'   for the `"phi"`-family excess-variance estimators).
#' @param do_estimate_lfdr Whether to additionally estimate a
#'   kernel-smoothed local false discovery rate via [estimate_lfdr()].
#' @param lfdr_pi0_method,lfdr_pi0_fixed,lfdr_pi0_lambda,lfdr_M,lfdr_K,
#'   lfdr_kernel,lfdr_isotonic Local-fdr estimation settings, passed to
#'   [estimate_lfdr()].
#' @param min_var,min_C,eps Numerical floors.
#' @param verbose Currently unused (reserved).
#'
#' @return A list of edge-level vectors and edge x donor matrices --
#'   see the file-level output field glossary above.
#' @keywords internal
fit_EWVDS_decomposition <- function(
  Y_mat,
  S2_mat,
  donor_ncells,
  noise_method = c("fisher_z", "bootstrap"),
  bootstrap_weight_method = c("fisher_z", "bootstrap_iv"),
  compute_Q0_sim = FALSE,
  null_method = c("pseudo_donor", "gaussian"),
  pseudo_donor_mode = c("permutation", "resampling"),
  Z_cell_mat = NULL,
  Z_g_edge = NULL,
  idx = NULL,
  pseudo_donor_replace = FALSE,
  use_pairwise = FALSE,
  n_cores = 1,
  n_null_sim = 999,
  null_seed = 2,
  evidence_p_method = c("Q0_DL", "Q0_sim"),
  excess_variance_method = c(
    "DL_tau",
    "DL_phi",
    "sim_tau",
    "sim_phi"
  ),
  phi_alpha = 0.5,
  do_estimate_lfdr = TRUE,
  lfdr_pi0_method = c("storey", "fixed"),
  lfdr_pi0_fixed = NULL,
  lfdr_pi0_lambda = 0.5,
  lfdr_M = 512,
  lfdr_K = 31,
  lfdr_kernel = c("triangular","uniform"),
  lfdr_isotonic = c("pava","cummax","none"),

  min_var = 1e-10,
  min_C = 1e-12,
  eps = 1e-12,

  verbose = TRUE
){
  lfdr_kernel <- match.arg(lfdr_kernel)
  lfdr_isotonic <- match.arg(lfdr_isotonic)
  noise_method <- match.arg(noise_method)
  bootstrap_weight_method <- match.arg(bootstrap_weight_method)
  null_method <- match.arg(null_method)
  pseudo_donor_mode <- match.arg(pseudo_donor_mode)
  evidence_p_method <- match.arg(evidence_p_method)
  excess_variance_method <- match.arg(excess_variance_method)
  lfdr_pi0_method <- match.arg(lfdr_pi0_method)
  lfdr_kernel <- match.arg(lfdr_kernel)

  if (!compute_Q0_sim && evidence_p_method == "Q0_sim") {
    stop(
      "evidence_p_method = \"Q0_sim\" requires compute_Q0_sim = TRUE ",
      "(the simulation-based null is never computed when compute_Q0_sim ",
      "is FALSE, which would otherwise leave p_Q_e silently NA). Either ",
      "set compute_Q0_sim = TRUE, or set evidence_p_method = \"Q0_DL\"."
    )
  }

  Y_mat <- as.matrix(Y_mat)
  S2_mat <- as.matrix(S2_mat)

  if(is.null(colnames(Y_mat))) stop("Y_mat must have donor colnames.")
  if(is.null(colnames(S2_mat))) colnames(S2_mat) <- colnames(Y_mat)
  if(!identical(colnames(S2_mat), colnames(Y_mat))) stop("colnames(S2_mat) must match colnames(Y_mat).")

  if(is.null(names(donor_ncells))){
  if(length(donor_ncells) != ncol(Y_mat)) stop("Unnamed donor_ncells length does not match ncol(Y_mat).")
    names(donor_ncells) <- colnames(Y_mat)
  }

  donor_ncells <- donor_ncells[colnames(Y_mat)]
  if(any(is.na(donor_ncells))) stop("donor_ncells is missing one or more donors in colnames(Y_mat).")

  S2_mat <- pmax(S2_mat, min_var)
  colnames(S2_mat) <- colnames(Y_mat)

  E <- nrow(Y_mat)
  D <- ncol(Y_mat)

  if(noise_method == "fisher_z" || bootstrap_weight_method == "fisher_z"){
    w_d <- as.numeric(donor_ncells[colnames(Y_mat)] - 3)
    W0 <- matrix(rep(w_d, each = E), nrow = E, ncol = D)
    colnames(W0) <- colnames(Y_mat)
  } else {
    W0 <- 1 / S2_mat
    colnames(W0) <- colnames(Y_mat)
  }

  W0 <- pmax(W0, eps)
  sumW0 <- rowSums(W0)
  A0 <- W0 / sumW0
  alpha0_e <- rowSums(W0 * Y_mat) / sumW0
  resid0 <- Y_mat - alpha0_e
  Q0_e <- rowSums(W0 * resid0^2)
  C0_e <- sumW0 - rowSums(W0^2) / sumW0
  C0_e <- pmax(C0_e, min_C)
  df_Q <- D - 1
  p_Q0_DL_e <- pchisq(Q0_e, df = df_Q, lower.tail = FALSE)
  tau2_DL_e <- pmax((Q0_e - df_Q) / C0_e, 0)
  r_DL_e <- Q0_e / df_Q - 1
  phi2_DL_e <- pmax(
    (Q0_e - qchisq(1 - phi_alpha, df = df_Q)) / C0_e,
    0
  )

  # =========================================================
  # Simulation branch (optional)
  # =========================================================

  p_Q0_sim_e <- rep(NA_real_, E)
  tau2_sim_e <- rep(NA_real_, E)
  r_sim_e <- rep(NA_real_, E)
  phi2_sim_e <- rep(NA_real_, E)
  Q0_null_mean_e <- rep(NA_real_, E)
  Q0_null_q_e <- rep(NA_real_, E)

  if(compute_Q0_sim){
    if(null_method == "gaussian"){
      set.seed(null_seed)
      Q0_null_mat <- matrix(
        0,
        nrow = E,
        ncol = n_null_sim
      )
      for(b in seq_len(n_null_sim)){
        Y_null <- matrix(
          rnorm(
            length(S2_mat),
            mean = 0,
            sd = sqrt(as.vector(S2_mat))
          ),
          nrow = E,
          ncol = D
        )
        alpha_null <- rowSums(W0 * Y_null) / sumW0
        resid_null <- Y_null - alpha_null
        Q0_null_mat[, b] <- rowSums(
          W0 * resid_null^2
        )
      }
    } else if(null_method == "pseudo_donor"){
      if(is.null(Z_cell_mat)) stop("Z_cell_mat required for pseudo_donor null.")
      if(is.null(Z_g_edge)) stop("Z_g_edge required for pseudo_donor null.")
      if(is.null(idx)) stop("idx required for pseudo_donor null.")
      pseudo_res <- simulate_Q0_null_pseudodonor(
        Z_cell_mat = Z_cell_mat,
        Z_g_edge = Z_g_edge,
        donor_ncells = donor_ncells,
        idx = idx,
        W0 = W0,
        n_sim = n_null_sim,
        mode = pseudo_donor_mode,
        replace = pseudo_donor_replace,
        seed = null_seed,
        use_pairwise = use_pairwise,
        eps = eps,
        n_cores_sim = n_cores
      )
      Q0_null_mat <- pseudo_res$Q0_null_mat
    } else {
      stop("Unknown null_method.")
    }
    Q0_null_mean_e <- rowMeans(Q0_null_mat, na.rm = TRUE)
    Q0_null_q_e <- apply(
      Q0_null_mat,
      1,
      quantile,
      probs = 1 - phi_alpha,
      na.rm = TRUE
    )
    p_Q0_sim_e <- (
      1 +
      rowSums(Q0_null_mat >= Q0_e, na.rm = TRUE)
    ) / (
      n_null_sim + 1
    )

    tau2_sim_e <- pmax(
      (
        Q0_e -
        Q0_null_mean_e
      ) / C0_e,
      0
    )
    r_sim_e <- Q0_e / Q0_null_mean_e - 1

    phi2_sim_e <- pmax(
      (
        Q0_e -
        Q0_null_q_e
      ) / C0_e,
      0
    )

    rm(Q0_null_mat)
    gc()
  }

  # ---------------------------------------------------------
  # choose evidence
  # ---------------------------------------------------------

  p_evidence_e <- if(
    evidence_p_method == "Q0_DL"
  ){
    p_Q0_DL_e
  } else {
    p_Q0_sim_e
  }

  q_Q_BH_e <- rep(NA_real_, E)
  q_Q_storey_e <- rep(NA_real_, E)
  pi0_Q_storey <- NA_real_

  ok_q <- is.finite(p_evidence_e) & !is.na(p_evidence_e) & p_evidence_e >= 0 & p_evidence_e <= 1

  if(sum(ok_q) > 0){
    q_Q_BH_e[ok_q] <- p.adjust(p_evidence_e[ok_q], method = "BH")
    qfit <- qvalue::qvalue(p = p_evidence_e[ok_q], lambda = lfdr_pi0_lambda)
    pi0_Q_storey <- qfit$pi0
    q_Q_storey_e[ok_q] <- qfit$qvalues
  }

  # ---------------------------------------------------------
  # lfdr
  # ---------------------------------------------------------

  if(do_estimate_lfdr){
    lfdr_fit <- estimate_lfdr(
      p = p_evidence_e,
      evidence_p_method = evidence_p_method,
      n_sim = n_null_sim,
      M = lfdr_M,
      K = lfdr_K,
      kernel = lfdr_kernel,
      pi0_method = lfdr_pi0_method,
      pi0_fixed = lfdr_pi0_fixed,
      pi0_lambda = lfdr_pi0_lambda,
      isotonic = lfdr_isotonic,
      eps = eps
    )
    lfdr_e <- lfdr_fit$lfdr
    rho_e <- lfdr_fit$rho_e
  } else {
    lfdr_fit <- NULL
    lfdr_e <- rep(NA_real_, E)
    rho_e <- 1 - p_evidence_e
  }
  rho_e <- pmin(pmax(rho_e, 0), 1)

  # ---------------------------------------------------------
  # choose excess variance
  # ---------------------------------------------------------

  phi2_e <- switch(
    excess_variance_method,
    "DL_phi" = phi2_DL_e,
    "sim_phi" = phi2_sim_e,
    "DL_tau" = tau2_DL_e,
    "sim_tau" = tau2_sim_e
  )

  tau2_e <- switch(
    excess_variance_method,
    "DL_phi" = tau2_DL_e,
    "sim_phi" = tau2_sim_e,
    "DL_tau" = tau2_DL_e,
    "sim_tau" = tau2_sim_e
  )

  phi2_e <- pmax(phi2_e, 0)

  # ---------------------------------------------------------
  # Step 3: one-step geometry
  # ---------------------------------------------------------

  W1 <- 1 / pmax(S2_mat + tau2_e, eps)
  sumW1 <- rowSums(W1)
  alpha1_e <- rowSums(W1 * Y_mat) / sumW1
  resid1 <- Y_mat - alpha1_e
  Q1_e <- rowSums(W1 * resid1^2)
  C1_e <- sumW1 - rowSums(W1^2) / sumW1

  # ---------------------------------------------------------
  # alpha test
  # ---------------------------------------------------------

  se_alpha_e <- sqrt(1 / sumW1)
  z_alpha_e <- alpha1_e / pmax(se_alpha_e, eps)
  p_alpha_e <- 2 * pnorm(-abs(z_alpha_e))

  # ---------------------------------------------------------
  # shrinkage
  # ---------------------------------------------------------

  # NOTE: reliability_mat (and hence lambda_DL_mat below) rises toward 1
  # as a donor's own sampling variance S2_mat falls -- i.e. as that
  # donor's cell count grows. This means Y_shrunk/Y_DL for a HIGH-cell-
  # count donor leans more on that donor's own raw per-edge deviation
  # (resid1), while a LOW-cell-count donor gets pulled harder toward the
  # pooled, many-donor average alpha1_e. One consequence is that the
  # run-to-run (randomized-PIT jitter) noise of the FINAL SHRUNK estimate
  # (R_DL_mat) can *increase* with donor cell count -- the opposite
  # direction from the RAW per-donor correlation's own noise (which
  # decreases with cell count, as ordinary sampling theory would
  # predict). This is expected: a high-reliability donor's shrunk
  # estimate is allowed to track its own (noisier, less-shrunk) raw
  # value more closely, while a low-reliability donor's own jitter noise
  # gets proportionally suppressed by the small lambda it's assigned.
  # The absolute scale stays far below the raw correlation's own noise
  # throughout, because alpha1_e itself pools ~D donors' worth of
  # independent jitter, and lambda is well below 1 for any edge whose
  # true heterogeneity phi2_e is small relative to S2_mat (true for most
  # edges; only edges with large real heterogeneity push lambda close to
  # 1 even for the best-sampled donors).
  reliability_mat <- phi2_e / pmax(phi2_e + S2_mat, eps)
  lambda_mat      <- reliability_mat * rho_e
  u_hat_mat       <- lambda_mat * resid1
  Y_shrunk        <- alpha1_e + u_hat_mat

  # DL versions (no rho_e weighting)
  lambda_DL_mat   <- reliability_mat             # lambda^DL = phi^2/(phi^2+sigma^2)
  u_hat_DL_mat    <- lambda_DL_mat * resid1      # u_hat^DL = lambda^DL * (y - alpha_hat)
  Y_DL            <- alpha1_e + u_hat_DL_mat     # y^DL = alpha_hat + u_hat^DL


  list(
    phi2_e = phi2_e,
    tau2_e = tau2_e,
    rho_e = rho_e,
    lfdr_e = lfdr_e,
    W0 = W0,
    W1 = W1,

    Q0_e = Q0_e,
    C0_e = C0_e,
    p_Q_e = p_evidence_e,
    q_Q_BH_e = q_Q_BH_e,
    q_Q_storey_e = q_Q_storey_e,
    pi0_Q_storey = pi0_Q_storey,
    p_Q0_DL_e = p_Q0_DL_e,

    tau2_DL_e = tau2_DL_e,
    r_DL_e = r_DL_e,
    phi2_DL_e = phi2_DL_e,

    p_Q0_sim_e = p_Q0_sim_e,
    tau2_sim_e = tau2_sim_e,
    r_sim_e = r_sim_e,
    phi2_sim_e = phi2_sim_e,
    Q0_null_mean_e = Q0_null_mean_e,
    Q0_null_q_e = Q0_null_q_e,

    alpha0_e = alpha0_e,
    alpha1_e = alpha1_e,
    alpha_e = alpha1_e,

    Q1_e = Q1_e,
    C1_e = C1_e,

    se_alpha_e = se_alpha_e,
    z_alpha_e = z_alpha_e,
    p_alpha_e = p_alpha_e,

    lambda_mat = lambda_mat,
    reliability_mat = reliability_mat,
    u_hat_mat = u_hat_mat,

    lambda_DL_mat = lambda_DL_mat,
    u_hat_DL_mat = u_hat_DL_mat,

    Y_mat = Y_mat,
    Y_shrunk = Y_shrunk,
    Y_DL = Y_DL,
    S2_mat = S2_mat,

    lfdr_fit = lfdr_fit,
    f_p_e = if(do_estimate_lfdr) lfdr_fit$f_p else rep(NA_real_, E)
  )
}

#' Remove leading principal components from the donor-deviation matrix
#'
#' Removes the leading `n_pc` principal components (computed across
#' donors) from the edge x donor deviation matrix, to strip out
#' low-rank structure shared by many edges before the heterogeneity
#' decomposition. Pure: does not write files or modify its input
#' in-place.
#'
#' This is an exploratory option, reached only via `do_pc_removal =
#' TRUE` in [run_EWVDS_pipeline()]/[run_ewvds_stage()]. It is not
#' described in the associated paper, and no analysis there sets
#' `do_pc_removal = TRUE` -- treat it as unverified if you use it.
#'
#' @param Y_mat Edge x donor deviation matrix.
#' @param n_pc Number of leading principal components to remove.
#' @param center Whether to center `Y_mat` (across donors) before PCA.
#'
#' @return `Y_mat` with the leading `n_pc` components removed, same
#'   dimnames, with an added `"pc_removed"` attribute recording how many
#'   components were actually removed (`<= n_pc`, bounded by the number
#'   of donors/edges available).
#' @keywords internal
remove_EWVDS_pcs <- function(
  Y_mat,
  n_pc = 5,
  center = TRUE
) {
  if (!is.matrix(Y_mat)) Y_mat <- as.matrix(Y_mat)
  if (any(!is.finite(Y_mat))) {
    stop("Y_mat contains NA/NaN/Inf before PC removal.")
  }
  if (n_pc < 1) {
    stop("n_pc must be >= 1 when do_pc_removal = TRUE.")
  }

  Y_centered <- scale(Y_mat, center = center, scale = FALSE)

  pc <- prcomp(Y_centered, center = FALSE, scale. = FALSE)
  k <- min(as.integer(n_pc), ncol(pc$x), ncol(pc$rotation))

  if (k < 1) {
    stop("No valid PC is available for PC removal.")
  }

  pcs <- pc$x[, seq_len(k), drop = FALSE]
  loadings <- pc$rotation[, seq_len(k), drop = FALSE]
  Y_pc <- pcs %*% t(loadings)

  Y_resid <- Y_centered - Y_pc
  rownames(Y_resid) <- rownames(Y_mat)
  colnames(Y_resid) <- colnames(Y_mat)

  attr(Y_resid, "pc_removed") <- k
  Y_resid
}

#' Bootstrap estimate of per-edge, within-donor sampling variance
#'
#' Estimates `sigma^2_de = Var_bootstrap(FisherZ(cor_de))` for every
#' edge `e` within one donor `d`, by resampling that donor's cells with
#' replacement. An alternative to the closed-form Fisher-z sampling
#' variance (`1 / (n_d - 3)`) used by [load_EWVDS_inputs()] when
#' `noise_method = "fisher_z"`.
#'
#' Reached only via `noise_method = "bootstrap"`, which the associated
#' paper's analyses do not use (they use `"fisher_z"` throughout) --
#' kept as a legitimate, documented method option, but not independently
#' tested or validated by any analysis in that paper.
#'
#' @param Z_donor Cells x genes matrix of this donor's Gaussian-copula
#'   scores.
#' @param idx Two-column edge index matrix (see [get_upper_idx()]).
#' @param n_boot Number of bootstrap replicates.
#' @param seed Random seed.
#' @param use_pairwise If `TRUE`, correlations use
#'   `"pairwise.complete.obs"` instead of `"everything"`.
#' @param eps Clip value before the Fisher-z transform.
#'
#' @return Numeric vector, one bootstrap variance estimate per edge.
#' @keywords internal
bootstrap_edge_variance <- function(
  Z_donor,
  idx,
  n_boot = 50,
  seed = 1,
  use_pairwise = FALSE,
  eps = 1e-6
) {
  set.seed(seed)

  n <- nrow(Z_donor)
  E <- nrow(idx)
  cor_use <- if (use_pairwise) "pairwise.complete.obs" else "everything"

  Z_boot_mat <- matrix(NA_real_, nrow = n_boot, ncol = E)

  for (b in seq_len(n_boot)) {
    samp <- sample.int(n, size = n, replace = TRUE)
    Sigma_b <- cor(Z_donor[samp, , drop = FALSE], use = cor_use)
    corr_b <- Sigma_b[cbind(idx[, 1], idx[, 2])]
    Z_boot_mat[b, ] <- fisher_z(corr_b, eps = eps)
  }

  apply(Z_boot_mat, 2, var, na.rm = TRUE)
}

# ============================================================
# Pipeline wrapper
# ============================================================

#' Load inputs, fit, and save Stage 2's complete output for one
#' configuration
#'
#' Wraps [load_EWVDS_inputs()], optionally [remove_EWVDS_pcs()],
#' [fit_EWVDS_decomposition()], and [save_EWVDS_donor_networks()]
#' into a single call, writing every intermediate and final result under
#' `COMB_PATH/<output_subdir>/` (named from `noise_method`,
#' `do_pc_removal`, and `pipeline_tag`).
#'
#' @inheritParams load_EWVDS_inputs
#' @inheritParams fit_EWVDS_decomposition
#' @param ensure_psd See [save_EWVDS_donor_networks()]. Default `FALSE`.
#' @param do_pc_removal,n_pc See [remove_EWVDS_pcs()]; exploratory,
#'   not used in the associated paper.
#' @param bootstrap_weight_method,n_boot,seed_bootstrap See
#'   [bootstrap_edge_variance()]; `noise_method = "bootstrap"` is not
#'   used in the associated paper.
#' @param lfdr_M,lfdr_K,lfdr_kernel,lfdr_pi0_method,lfdr_pi0_fixed,lfdr_pi0_lambda,lfdr_isotonic See [estimate_lfdr()].
#' @param min_var,min_C,eps Numerical floors used in Stage 2's
#'   random-effects decomposition (see
#'   [fit_EWVDS_decomposition()]); `eps` is also used as the
#'   clip value before Fisher-z transforms in [load_EWVDS_inputs()].
#' @param force_recompute_stage2 If `FALSE` (default), reuse cached
#'   `Sigma_donor.rds`/`bootstrap_var_fisherZ.rds`/`Z_randomized.rds`
#'   under `COMB_PATH` where present (see [load_EWVDS_inputs()]). This
#'   only affects Stage 2's own cache -- it has no effect on Stage 1a/1b.
#' @param pipeline_tag String appended to the output subdirectory name.
#' @param verbose Currently unused (reserved).
#'
#' @return Invisibly, a list with `inputs` (from [load_EWVDS_inputs()]),
#'   `fit` (from [fit_EWVDS_decomposition()]), and `base_dir`.
#'
#' @section Output:
#' Writes `EWVDS_inputs.rds`, `EWVDS_full_intermediate.rds`,
#' `EWVDS_result.rds`, and per-donor `Sigma_EWVDS.rds`/`Sigma_DL.rds`
#' under `COMB_PATH/<output_subdir>/`. See the file-level output field
#' glossary above for what each field means.
#' @keywords internal
run_EWVDS_pipeline <- function(
  COMB_PATH,
  donors,
  noise_method = c("fisher_z", "bootstrap"),
  do_pc_removal = FALSE,
  n_pc = 5,
  bootstrap_weight_method = c("fisher_z", "bootstrap_iv"),
  n_boot = 50,
  seed_bootstrap = 42,
  compute_Q0_sim = FALSE,
  null_method = c("pseudo_donor", "gaussian"),
  pseudo_donor_mode = c("permutation", "resampling"),
  pseudo_donor_replace = FALSE,
  n_cores = 1,
  n_null_sim = 999,
  seed_Q0_sim = 43,
  evidence_p_method = c("Q0_DL","Q0_sim"),
  excess_variance_method = c("DL_tau", "DL_phi", "sim_tau", "sim_phi"),
  phi_alpha = 0.50,
  do_estimate_lfdr = TRUE,
  lfdr_M = NULL,
  lfdr_K = 31,
  lfdr_kernel = c("triangular", "uniform"),
  lfdr_pi0_method = c("storey", "fixed"),
  lfdr_pi0_fixed = NULL,
  lfdr_pi0_lambda = 0.5,
  lfdr_isotonic = c("pava", "cummax", "none"),
  force_recompute_stage2 = FALSE,
  use_pairwise = FALSE,
  min_var = 1e-10,
  min_C = 1e-12,
  eps = 1e-12,
  pipeline_tag = "pipeline",
  verbose = TRUE,
  ensure_psd = FALSE
){

  noise_method <- match.arg(noise_method)
  bootstrap_weight_method <- match.arg(bootstrap_weight_method)
  null_method <- match.arg(null_method)
  pseudo_donor_mode <- match.arg(pseudo_donor_mode)
  evidence_p_method <- match.arg(evidence_p_method)
  excess_variance_method <- match.arg(excess_variance_method)
  lfdr_kernel <- match.arg(lfdr_kernel)
  lfdr_pi0_method <- match.arg(lfdr_pi0_method)
  lfdr_isotonic <- match.arg(lfdr_isotonic)

  if (is.null(lfdr_M)) lfdr_M <- min(n_null_sim, 1000L)

  output_subdir <- "donor_network_EWVDS"
  if(noise_method == "bootstrap"){
    output_subdir <- paste0(output_subdir, "_bootstrap")
    if(bootstrap_weight_method == "fisher_z") output_subdir <- paste0(output_subdir, "_IVfisher")
    if(bootstrap_weight_method == "bootstrap_iv") output_subdir <- paste0(output_subdir, "_IVbootstrap")
  }

  if(do_pc_removal) output_subdir <- paste0(output_subdir, "_PC", n_pc)
  if(!is.null(pipeline_tag)) output_subdir <- paste0(output_subdir, "_", pipeline_tag)

  inputs <- load_EWVDS_inputs(
    COMB_PATH=COMB_PATH,
    donors=donors,
    n_boot=n_boot,
    noise_method=noise_method,
    seed_bootstrap=seed_bootstrap,
    force_recompute_stage2=force_recompute_stage2,
    use_pairwise=use_pairwise,
    eps=eps,
    n_cores=n_cores,
    output_subdir=output_subdir
  )

  if(do_pc_removal) inputs$Y_mat <- remove_EWVDS_pcs(Y_mat=inputs$Y_mat, n_pc=n_pc, center=TRUE)

  fit <- fit_EWVDS_decomposition(
    Y_mat=inputs$Y_mat,
    S2_mat=inputs$S2_mat,
    donor_ncells=inputs$donor_ncells,
    noise_method=noise_method,
    bootstrap_weight_method=bootstrap_weight_method,
    compute_Q0_sim=compute_Q0_sim,
    null_method=null_method,
    pseudo_donor_mode=pseudo_donor_mode,
    Z_cell_mat=inputs$Z_cell_mat,
    Z_g_edge=inputs$Z_g_edge,
    idx=inputs$idx,
    pseudo_donor_replace=pseudo_donor_replace,
    use_pairwise=use_pairwise,
    n_cores=n_cores,
    n_null_sim=n_null_sim,
    null_seed=seed_Q0_sim,
    evidence_p_method=evidence_p_method,
    excess_variance_method=excess_variance_method,
    phi_alpha=phi_alpha,
    do_estimate_lfdr=do_estimate_lfdr,
    lfdr_M=lfdr_M,
    lfdr_K=lfdr_K,
    lfdr_kernel=lfdr_kernel,
    lfdr_pi0_method=lfdr_pi0_method,
    lfdr_pi0_fixed=lfdr_pi0_fixed,
    lfdr_pi0_lambda=lfdr_pi0_lambda,
    lfdr_isotonic=lfdr_isotonic,
    min_var=min_var,
    min_C=min_C,
    eps=eps,
    verbose=verbose
  )

  base_dir    <- inputs$base_dir
  # R_EWVDS_mat is derived from fit$Y_shrunk, which additionally applies
  # rho_e = 1-lfdr_e weighting on top of the reliability-based shrinkage
  # (see Y_shrunk vs Y_DL above). Downstream evaluation code should read
  # R_DL_mat (and lambda_DL_mat), not R_EWVDS_mat/lambda_mat, unless it
  # specifically wants the local-FDR-weighted variant -- see the
  # file-level output field glossary above.
  Z_EWVDS_mat <- sweep(fit$Y_shrunk, 1, inputs$Z_g_edge, "+")
  R_EWVDS_mat <- inv_fisher_z(Z_EWVDS_mat)
  Z_DL_mat    <- sweep(fit$Y_DL,     1, inputs$Z_g_edge, "+")
  R_DL_mat    <- inv_fisher_z(Z_DL_mat)

  edge_df <- data.frame(
    edge_id = seq_len(nrow(inputs$idx)),
    g1 = inputs$idx[, 1],
    g2 = inputs$idx[, 2],
    gene1 = inputs$gene_names[inputs$idx[, 1]],
    gene2 = inputs$gene_names[inputs$idx[, 2]]
  )
  edge_df$p_Q_e <- fit$p_Q_e
  edge_df$q_Q_BH_e <- fit$q_Q_BH_e
  edge_df$q_Q_storey_e <- fit$q_Q_storey_e
  edge_df$lfdr_e <- fit$lfdr_e
  edge_df$rho_e <- fit$rho_e
  edge_df$phi2_e <- fit$phi2_e
  edge_df$tau2_e <- fit$tau2_e

  saveRDS(c(inputs, fit, list(Z_EWVDS_mat=Z_EWVDS_mat, R_EWVDS_mat=R_EWVDS_mat,
                              Z_DL_mat=Z_DL_mat, R_DL_mat=R_DL_mat,
                              edge_df=edge_df)), file.path(base_dir, "EWVDS_full_intermediate.rds"))
  saveRDS(
    list(
      phi2_e=fit$phi2_e, tau2_e=fit$tau2_e, rho_e=fit$rho_e,
      lfdr_e=fit$lfdr_e, f_p_e=fit$f_p_e, lfdr_fit=fit$lfdr_fit,
      Q0_e=fit$Q0_e, C0_e=fit$C0_e, Q1_e=fit$Q1_e, C1_e=fit$C1_e,
      p_Q0_DL_e=fit$p_Q0_DL_e, tau2_DL_e=fit$tau2_DL_e, r_DL_e=fit$r_DL_e, phi2_DL_e=fit$phi2_DL_e,
      p_Q0_sim_e=fit$p_Q0_sim_e, tau2_sim_e=fit$tau2_sim_e, r_sim_e=fit$r_sim_e, phi2_sim_e=fit$phi2_sim_e,
      Q0_null_mean_e=fit$Q0_null_mean_e, Q0_null_q_e=fit$Q0_null_q_e,
      p_Q_e=fit$p_Q_e, q_Q_BH_e=fit$q_Q_BH_e, q_Q_storey_e=fit$q_Q_storey_e, pi0_Q_storey=fit$pi0_Q_storey,
      evidence_p_method=evidence_p_method, null_method=null_method, pseudo_donor_mode=pseudo_donor_mode, pseudo_donor_replace=pseudo_donor_replace, n_null_sim=n_null_sim, n_cores=n_cores,
      noise_method=noise_method, bootstrap_weight_method=bootstrap_weight_method, seed_Q0_sim=seed_Q0_sim, phi_alpha=phi_alpha, excess_variance_method=excess_variance_method,
      COMB_PATH=COMB_PATH, use_pairwise=use_pairwise, do_pc_removal=do_pc_removal, n_pc=n_pc,
      lfdr_M=lfdr_M, lfdr_K=lfdr_K, lfdr_kernel=lfdr_kernel, lfdr_pi0_method=lfdr_pi0_method, lfdr_pi0_fixed=lfdr_pi0_fixed, lfdr_pi0_lambda=lfdr_pi0_lambda, lfdr_isotonic=lfdr_isotonic,
      alpha0_e=fit$alpha0_e, alpha1_e=fit$alpha1_e, alpha_e=fit$alpha_e,
      se_alpha_e=fit$se_alpha_e, z_alpha_e=fit$z_alpha_e, p_alpha_e=fit$p_alpha_e,
      S2_mat=fit$S2_mat, W0=fit$W0, W1=fit$W1, lambda_mat=fit$lambda_mat, reliability_mat=fit$reliability_mat,
      lambda_DL_mat=fit$lambda_DL_mat, u_hat_DL_mat=fit$u_hat_DL_mat,
      Y_mat=fit$Y_mat, Y_shrunk=fit$Y_shrunk, Y_DL=fit$Y_DL,
      Z_EWVDS_mat=Z_EWVDS_mat, R_EWVDS_mat=R_EWVDS_mat, u_hat_mat=fit$u_hat_mat,
      Z_DL_mat=Z_DL_mat, R_DL_mat=R_DL_mat,
      idx=inputs$idx, donors=inputs$donors, donor_ncells=inputs$donor_ncells, n_null_sim=n_null_sim,
      pipeline_version=pipeline_tag, edge_df=edge_df, gene_names=inputs$gene_names
    ),
    file.path(base_dir, "EWVDS_result.rds")
  )

  save_EWVDS_donor_networks(base_dir=base_dir, donors=inputs$donors, Sigma_global_Z=inputs$Sigma_global_Z, Z_g_edge=inputs$Z_g_edge, Y_shrunk=fit$Y_shrunk, idx=inputs$idx, gene_names=inputs$gene_names, Y_DL=fit$Y_DL, ensure_psd=ensure_psd)

  invisible(list(inputs=inputs, fit=fit, base_dir=base_dir))
}

#' Run Stage 2 (donor heterogeneity estimation and shrinkage)
#'
#' Thin wrapper around [run_EWVDS_pipeline()]: auto-detects `donors`
#' from `comb_path/donors/` if not supplied, then runs the full Stage 2
#' decomposition and shrinkage. See [run_EWVDS_pipeline()] and
#' [fit_EWVDS_decomposition()] for parameter details.
#'
#' @inheritParams run_EWVDS_pipeline
#' @param comb_path Root output directory, same as passed to
#'   [run_marginal_stage()] / [run_donor_split()].
#' @param bootstrap_weight_method See [bootstrap_edge_variance()];
#'   `noise_method = "bootstrap"` is not used in the associated paper.
#' @param seed_bootstrap Random seed for the bootstrap variance path.
#' @param seed_Q0_sim Random seed for the pseudo-donor null simulation.
#' @param lfdr_pi0_method,lfdr_pi0_lambda,lfdr_K,lfdr_kernel See
#'   [estimate_lfdr()].
#' @param min_var,min_C,eps Numerical floors used in Stage 2's
#'   random-effects decomposition.
#'
#' @return Invisibly, the result of [run_EWVDS_pipeline()].
#' @export
run_ewvds_stage <- function(
  comb_path,
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
  eps                     = 1e-12,
  pipeline_tag            = "pipeline_DL_noisefisher",
  n_cores                 = 30L,
  force_recompute_stage2         = FALSE,
  ensure_psd              = FALSE
) {
  if (is.null(donors)) {
    donors <- list.dirs(file.path(comb_path, "donors"),
                        full.names = FALSE, recursive = FALSE)
    message(sprintf("[Stage 2: heterogeneity + shrinkage] Auto-detected %d donors from %s/donors/",
                    length(donors), comb_path))
  }

  message(sprintf("[Stage 2: heterogeneity + shrinkage] Running on %d donors ...", length(donors)))

  res <- run_EWVDS_pipeline(
    COMB_PATH               = comb_path,
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
    eps                     = eps,
    pipeline_tag            = pipeline_tag,
    n_cores                 = n_cores,
    force_recompute_stage2         = force_recompute_stage2,
    ensure_psd              = ensure_psd
  )

  message("[Stage 2: heterogeneity + shrinkage] Done.")
  invisible(res)
}
