# -----------------------------
# Pseudo-donor null utilities, for the simulation-based (as opposed to
# asymptotic chi-squared) null distribution of the heterogeneity
# statistic Q_e.
#
# Reached only via compute_Q0_sim = TRUE + null_method = "pseudo_donor":
# not used by the associated paper's main reported results (which use
# the asymptotic DL null throughout), but used by its DL-vs-simulation
# robustness check (see the appendix section comparing the two nulls).
#
# Mechanism note: all B replicates reuse the SAME observed Z_cell_mat
# (one fixed randomized-PIT jitter realization) -- only cell-to-
# pseudo-donor-bucket membership is reshuffled per replicate (via
# sample.int(), mode "resampling" or "permutation"). This is NOT the
# same thing as redrawing the randomized-PIT jitter itself.
# -----------------------------

#' Assign cells to pseudo-donor buckets
#'
#' @param donor_ncells Named integer vector, one cell count per
#'   (pseudo-)donor.
#' @param n_cell_total Total number of cells to draw bucket membership
#'   from.
#' @param mode `"permutation"` (default, used in the associated paper:
#'   a single permutation of `sum(donor_ncells)` cells is split into
#'   consecutive, non-overlapping buckets; requires
#'   `sum(donor_ncells) <= n_cell_total`) or `"resampling"` (each
#'   pseudo-donor's cells are drawn independently, with or without
#'   replacement depending on `replace`, so the same cell can appear in
#'   more than one pseudo-donor; this biases `mean(Q_e^null)` slightly
#'   below its `D - 1` target).
#' @param replace Only used when `mode = "resampling"`.
#'
#' @return A named list of integer index vectors, one per entry of
#'   `donor_ncells`.
#' @keywords internal
generate_pseudodonor_indices <- function(
  donor_ncells,
  n_cell_total,
  mode = c("permutation", "resampling"),
  replace = FALSE
){
  mode <- match.arg(mode)

  if(mode == "resampling"){
    idx_list <- lapply(
      donor_ncells,
      function(n_d){
        sample.int(
          n = n_cell_total,
          size = n_d,
          replace = replace
        )
      }
    )
    names(idx_list) <- names(donor_ncells)
    return(idx_list)
  }

  if(mode == "permutation"){
    total_needed <- sum(donor_ncells)
    if(total_needed > n_cell_total){
      stop("Permutation mode requires sum(donor_ncells) <= n_cell_total.")
    }

    perm <- sample.int(
      n = n_cell_total,
      size = total_needed,
      replace = FALSE
    )

    idx_list <- split(
      perm,
      rep(seq_along(donor_ncells), donor_ncells)
    )

    names(idx_list) <- names(donor_ncells)
    return(idx_list)
  }
}

#' Edge-level deviations for one pseudo-donor assignment
#'
#' Given a pseudo-donor-to-cell-index assignment (from
#' [generate_pseudodonor_indices()]), computes each pseudo-donor's own
#' correlation matrix from its assigned cells, Fisher-z transforms it,
#' and subtracts the observed global correlation to get the same
#' edge x donor deviation structure Stage 2 operates on.
#'
#' @param Z_cell_mat Cells x genes matrix of Gaussian-copula scores (the
#'   single, fixed, observed jitter realization -- not redrawn per
#'   pseudo-donor or per replicate).
#' @param Z_g_edge Edge-level global correlation, Fisher-z scale.
#' @param pseudo_idx_list Named list of cell-index vectors, one per
#'   pseudo-donor (from [generate_pseudodonor_indices()]).
#' @param idx Two-column edge index matrix (see [get_upper_idx()]).
#' @param use_pairwise If `TRUE`, correlations use
#'   `"pairwise.complete.obs"` instead of `"everything"`.
#' @param eps Clip value before the Fisher-z transform.
#'
#' @return A list with `Z_null_mat` (edge x pseudo-donor, Fisher-z
#'   scale), `Z_g_edge`, and `Y_null_mat` (`Z_null_mat` minus
#'   `Z_g_edge`).
#' @keywords internal
compute_pseudodonor_Zmat <- function(
  Z_cell_mat,
  Z_g_edge,
  pseudo_idx_list,
  idx,
  use_pairwise = FALSE,
  eps = 1e-6
){
  cor_use <- if(use_pairwise) "pairwise.complete.obs" else "everything"
  D <- length(pseudo_idx_list)
  E <- nrow(idx)

  Z_null_mat <- matrix(
    NA_real_,
    nrow = E,
    ncol = D
  )

  for(d in seq_len(D)){

    Z_sub <- Z_cell_mat[
      pseudo_idx_list[[d]],
      ,
      drop = FALSE
    ]

    Sigma_d <- suppressWarnings(
      cor(
        Z_sub,
        use = cor_use
      )
    )

    Sigma_d[is.na(Sigma_d)] <- 0
    diag(Sigma_d) <- 1

    Sigma_d <- clip_cor(
      Sigma_d,
      eps = eps
    )

    Sigma_d_Z <- fisher_z(
      Sigma_d,
      eps = eps
    )

    Z_null_mat[, d] <- Sigma_d_Z[
      cbind(idx[,1], idx[,2])
    ]
  }

  colnames(Z_null_mat) <- names(pseudo_idx_list)
  if(length(Z_g_edge) != nrow(Z_null_mat)){
    stop("length(Z_g_edge) must match nrow(Z_null_mat).")
  }

  Y_null_mat <- sweep(
    Z_null_mat,
    1,
    Z_g_edge,
    "-"
  )

  list(
    Z_null_mat = Z_null_mat,
    Z_g_edge = Z_g_edge,
    Y_null_mat = Y_null_mat
  )
}

#' Simulation-based (pseudo-donor) null distribution of Q_e
#'
#' Draws `n_sim` independent pseudo-donor bucket assignments (reusing
#' the single observed `Z_cell_mat` jitter realization every time -- see
#' the file-level note above), and for each one recomputes the
#' heterogeneity statistic `Q_e` under the null hypothesis of no true
#' donor heterogeneity, using the same weights `W0` as the observed fit.
#'
#' @param Z_cell_mat Cells x genes matrix of Gaussian-copula scores.
#' @param Z_g_edge Edge-level global correlation, Fisher-z scale.
#' @param donor_ncells Named integer vector of (pseudo-)donor cell
#'   counts, matching `colnames(W0)`.
#' @param idx Two-column edge index matrix (see [get_upper_idx()]).
#' @param W0 Edge x donor weight matrix from the observed fit (see
#'   [fit_EWVDS_decomposition()]).
#' @param n_sim Number of null replicates.
#' @param mode,replace Passed to [generate_pseudodonor_indices()].
#' @param seed Random seed.
#' @param use_pairwise,eps Passed to [compute_pseudodonor_Zmat()].
#' @param n_cores_sim Number of parallel workers across replicates. With
#'   `n_cores_sim > 1`, replicates are drawn from L'Ecuyer-CMRG streams so
#'   the result is reproducible for a given `seed` and `n_cores_sim` (it
#'   depends on `n_cores_sim`, because replicates are prescheduled across
#'   workers). `n_cores_sim = 1` uses the caller's RNG kind, unchanged from
#'   earlier versions.
#'
#' @return A list with `Q0_null_mat`, an edge x `n_sim` matrix of
#'   null-replicate `Q_e` values.
#' @keywords internal
simulate_Q0_null_pseudodonor <- function(
  Z_cell_mat,
  Z_g_edge,
  donor_ncells,
  idx,
  W0,
  n_sim = 1000,
  mode = c("permutation", "resampling"),
  replace = FALSE,
  seed = 1,
  use_pairwise = FALSE,
  eps = 1e-6,
  n_cores_sim = 1
){

  mode <- match.arg(mode)
  if (n_cores_sim > 1) {
    # Under the default Mersenne-Twister, mclapply() re-seeds each forked
    # child from time and PID, so set.seed() alone left the parallel path
    # non-reproducible. L'Ecuyer-CMRG plus mc.reset.stream() ties the
    # children's streams to `seed`. The serial path is left as it was so
    # its output does not change.
    old_rngkind <- RNGkind()[1]
    on.exit(RNGkind(old_rngkind), add = TRUE)
    RNGkind("L'Ecuyer-CMRG")
  }
  set.seed(seed)
  if (n_cores_sim > 1) parallel::mc.reset.stream()

  E <- nrow(idx)
  D <- length(donor_ncells)

  if(is.null(colnames(W0))){
    stop("W0 must have donor colnames.")
  }

  if(is.null(names(donor_ncells))){
    stop("donor_ncells must have donor names.")
  }

  donor_ncells <- donor_ncells[colnames(W0)]

  if(any(is.na(donor_ncells))){
    stop("donor_ncells is missing one or more donors in colnames(W0).")
  }

  if(ncol(W0) != length(donor_ncells)){
    stop("ncol(W0) does not match length(donor_ncells).")
  }

  Q0_null_mat <- matrix(
    0,
    nrow = E,
    ncol = n_sim
  )

  sumW0 <- rowSums(W0)

  sim_list <- parallel::mclapply(
    seq_len(n_sim),
    function(b){

      pseudo_idx_list <- generate_pseudodonor_indices(
        donor_ncells = donor_ncells,
        n_cell_total = nrow(Z_cell_mat),
        mode = mode,
        replace = replace
      )

      names(pseudo_idx_list) <- names(donor_ncells)

      pseudo_res <- compute_pseudodonor_Zmat(
        Z_cell_mat = Z_cell_mat,
        Z_g_edge = Z_g_edge,
        pseudo_idx_list = pseudo_idx_list,
        idx = idx,
        use_pairwise = use_pairwise,
        eps = eps
      )

      Y_null <- pseudo_res$Y_null_mat

      if(is.null(colnames(Y_null))){
        stop("Y_null has no colnames; cannot align pseudo-donor null with W0.")
      }

      if(is.null(colnames(W0))){
        stop("W0 has no colnames; cannot align pseudo-donor null with Y_null.")
      }

      if(!setequal(colnames(Y_null), colnames(W0))){
        stop("Y_null and W0 have different donor names.")
      }

      Y_null <- Y_null[, colnames(W0), drop = FALSE]

      if(!identical(colnames(Y_null), colnames(W0))){
        stop("Failed to align Y_null columns to W0 columns.")
      }

      alpha0_null <- rowSums(
        W0 * Y_null
      ) / sumW0

      resid0_null <- Y_null - alpha0_null

      rowSums(
        W0 * resid0_null^2
      )
    },
    mc.cores = n_cores_sim
  )

  failed <- vapply(sim_list, inherits, logical(1), what = "try-error")
  if (any(failed)) {
    n_failed <- sum(failed)
    first_err <- conditionMessage(attr(sim_list[[which(failed)[1]]], "condition"))
    stop(sprintf(
      "simulate_Q0_null_pseudodonor: %d / %d simulation(s) failed.\nFirst error: %s",
      n_failed, n_sim, first_err
    ))
  }

  Q0_null_mat <- do.call(
    cbind,
    sim_list
  )

  list(
    Q0_null_mat = Q0_null_mat
  )
}
