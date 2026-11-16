# popcorn 0.1.1

First public release, accompanying the associated paper. The paper's
analysis code (<https://github.com/keleslab/popcorn-paper>) uses this
version.

## New features

* New `seed` argument (default `NULL`) in `run_coexp_pipeline()` and
  `run_marginal_stage()`. When given, Stage 1 runs `set.seed(seed)` and
  `parallel::mc.reset.stream()` under L'Ecuyer-CMRG, so the randomized
  PIT, and therefore the whole pipeline, is bit-reproducible, including
  across repeated calls in one R session and across `n_cores`. With the
  default `seed = NULL`, behavior is unchanged.
* New `ensure_psd` option (default `FALSE`) in `run_coexp_pipeline()`,
  `run_ewvds_stage()`, and `save_EWVDS_donor_networks()`. POPCORN
  shrinks each edge separately, so a donor's assembled correlation
  matrix `Sigma_DL.rds` is not guaranteed to be positive semi-definite.
  With `ensure_psd = TRUE`, each donor also gets `Sigma_DL_nearPD.rds`
  (the nearest correlation matrix, Higham 2002, via
  `Matrix::nearPD(corr = TRUE)`, applied only when the smallest
  eigenvalue is below `-psd_tol`; otherwise an exact copy), and Stage 2
  writes a per-donor summary `psd_check_DL.csv`. `Sigma_DL.rds` and
  every edge-level output (`EWVDS_result.rds`) are unchanged either way.
* New exported helper `ensure_psd_correlation()` for a single matrix.
* `celltype` no longer requires a real `colData(sce)` column. If the
  named column (default `"cell_type"`) does not exist, one is added with
  a single constant level. The argument only satisfies
  `scDesign3::construct_data()` and never subsets cells.

## Changed behavior and defaults

* The global correlation is now computed from the same Gaussian-copula
  scores as the per-donor correlations (PIT values clipped to
  `[1e-6, 1 - 1e-6]`). The `eps` argument of `run_donor_split()` and
  `eps_pit` of `run_coexp_pipeline()` are deprecated and ignored, with a
  warning. This shifts each edge's global correlation by a small,
  edge-specific constant; Stage 2's heterogeneity estimates and shrunk
  per-donor estimates do not depend on it.
* `donor_col` no longer has a default in `run_marginal_stage()`,
  `run_coexp_pipeline()`, and `run_donor_split()`; which column
  identifies donors must be given explicitly.
* `gene_names` defaults to `NULL` (all genes in `rownames(sce)`);
  `covariates` and `mu_formula` default to `"log1p_total_counts"`;
  `sigma_formula` defaults to `"1"`.
* `pseudo_donor_mode` defaults to `"permutation"` (disjoint
  pseudo-donors, as in the associated paper). This only matters when
  `compute_Q0_sim = TRUE`.
* `run_coexp_pipeline()`'s `pipeline_tag` defaults to
  `"pipeline_DL_noisefisher"`, matching `run_ewvds_stage()`.
* The pseudo-donor null simulation is reproducible under `seed_Q0_sim`
  when it runs on more than one core (L'Ecuyer-CMRG streams); its output
  then depends on the number of cores. Single-core output is unchanged.
* `Matrix` is now an imported package.

# popcorn 0.1.0

* Initial version: Stage 1 marginal fit (`run_marginal_stage()`),
  raw global and per-donor correlations (`run_donor_split()`), Stage 2
  random-effects decomposition and shrinkage (`run_ewvds_stage()`), and
  the end-to-end wrapper `run_coexp_pipeline()`.
* The asymptotic chi-squared null is the default for the heterogeneity
  test (`compute_Q0_sim = FALSE`, `evidence_p_method = "Q0_DL"`).
* `run_marginal_stage()` errors if `comb_path` already holds Stage 1
  output for a different `sce`; `evidence_p_method = "Q0_sim"` requires
  `compute_Q0_sim = TRUE`.
