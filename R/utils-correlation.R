# -----------------------------
# Basic correlation/Fisher-z utilities shared across Stage 1 and Stage 2.
# -----------------------------

#' Clip a correlation away from the +-1 boundary
#'
#' @param r Numeric vector of correlations.
#' @param eps Distance to keep `r` away from -1/1.
#'
#' @return Numeric vector, clipped to `[-1 + eps, 1 - eps]`.
#' @keywords internal
clip_cor <- function(r, eps = 1e-6) {
  pmin(pmax(r, -1 + eps), 1 - eps)
}

#' Fisher z-transform a correlation
#'
#' Clips `r` away from the +-1 boundary (see [clip_cor()]) before applying
#' `atanh()`, so this is safe to call on correlations saturated at +-1.
#'
#' @param r Numeric vector of correlations.
#' @param eps Distance to keep `r` away from -1/1 before transforming.
#'
#' @return Numeric vector of Fisher z-transformed values.
#' @export
fisher_z <- function(r, eps = 1e-6) {
  atanh(clip_cor(r, eps = eps))
}

#' Inverse Fisher z-transform
#'
#' @param z Numeric vector of Fisher z values.
#'
#' @return Numeric vector of correlations, `tanh(z)`.
#' @export
inv_fisher_z <- function(z) {
  tanh(z)
}

#' Upper-triangle edge indices of a square matrix
#'
#' @param Sigma A square matrix (only its dimensions are used).
#'
#' @return A two-column integer matrix of `(row, col)` indices into the
#'   upper triangle of `Sigma`, as returned by `which(..., arr.ind = TRUE)`.
#' @export
get_upper_idx <- function(Sigma) {
  which(upper.tri(Sigma), arr.ind = TRUE)
}

#' Rebuild a symmetric gene-by-gene matrix from edge-level values
#'
#' @param edge_values Numeric vector, one value per edge, in the same
#'   order as the rows of `idx`.
#' @param idx Two-column integer matrix of `(row, col)` indices, as
#'   returned by [get_upper_idx()].
#' @param gene_names Optional character vector of gene names to use as row
#'   and column names of the result.
#' @param diag_value Value to place on the diagonal (e.g. `1` for a
#'   correlation matrix).
#'
#' @return A symmetric `G x G` matrix, where `G = max(idx)`.
#' @export
edge_to_symmetric_mat <- function(edge_values, idx, gene_names = NULL, diag_value = 1) {
  G <- max(idx)
  M <- matrix(0, nrow = G, ncol = G)

  M[idx] <- edge_values
  M <- M + t(M)
  diag(M) <- diag_value

  if (!is.null(gene_names)) {
    rownames(M) <- gene_names
    colnames(M) <- gene_names
  }

  M
}

#' Project a correlation matrix to the nearest positive semi-definite one
#'
#' POPCORN estimates and shrinks each edge separately, so a donor-level
#' matrix assembled from its edge-level estimates is symmetric with unit
#' diagonal but is not guaranteed to be positive semi-definite (PSD).
#' This function checks the smallest eigenvalue and, only if it is below
#' `-tol`, replaces the matrix with the nearest correlation matrix in
#' Frobenius norm (Higham, 2002), computed by [Matrix::nearPD()] with
#' `corr = TRUE`. Matrices that are already PSD are returned unchanged.
#'
#' @param S Symmetric numeric matrix with unit diagonal.
#' @param tol Non-negative tolerance. `S` is treated as PSD when its
#'   smallest eigenvalue is at least `-tol`.
#'
#' @return A list with
#'   * `mat`: the PSD correlation matrix (`S` itself if no projection was
#'     needed), with `S`'s dimnames;
#'   * `projected`: whether a projection was applied;
#'   * `min_eig_before`, `min_eig_after`: smallest eigenvalue of `S` and
#'     of `mat`;
#'   * `max_abs_change`: largest absolute change over off-diagonal
#'     entries (`0` if not projected).
#'
#' @references Higham, N. J. (2002). Computing the nearest correlation
#'   matrix -- a problem from finance. *IMA Journal of Numerical
#'   Analysis*, 22(3), 329--343.
#' @export
ensure_psd_correlation <- function(S, tol = 1e-8) {
  stopifnot(is.matrix(S), nrow(S) == ncol(S), tol >= 0)
  min_eig <- function(M) min(eigen(M, symmetric = TRUE, only.values = TRUE)$values)

  ev0 <- min_eig(S)
  if (ev0 >= -tol) {
    return(list(mat = S, projected = FALSE, min_eig_before = ev0,
                min_eig_after = ev0, max_abs_change = 0))
  }

  np <- Matrix::nearPD(S, corr = TRUE)
  if (!isTRUE(np$converged)) {
    warning("Matrix::nearPD() did not converge; the returned matrix may not be fully projected.")
  }
  P <- as.matrix(np$mat)
  dimnames(P) <- dimnames(S)

  ut <- upper.tri(S)
  list(mat = P, projected = TRUE, min_eig_before = ev0,
       min_eig_after = min_eig(P), max_abs_change = max(abs(P[ut] - S[ut])))
}
