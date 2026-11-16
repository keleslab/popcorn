# ============================================================
# Stage 1: covariate-adjusted marginal model fit per gene, and the
# per-gene fitting loop that drives it (Methods, "Stage 1: Marginal
# Fitting and Raw Correlation").
# ============================================================

#' Check that required covariates are present in a SingleCellExperiment
#'
#' Validates that every name in `covariates` is a column of
#' `colData(sce)`, printing a short report and stopping with an
#' informative error if any are missing. Intended as a pre-flight check
#' before [run_marginal_stage()].
#'
#' @param sce A `SingleCellExperiment`.
#' @param covariates Character vector of `colData(sce)` column names to
#'   check for.
#'
#' @return Invisible `NULL`; called for its side effect (messages, or an
#'   error if a covariate is missing).
#' @keywords internal
check_covariates <- function(sce, covariates) {
  cd      <- SummarizedExperiment::colData(sce)
  missing <- setdiff(covariates, colnames(cd))
  cat("=== Checking covariates in colData(sce) ===\n")
  for (cov in covariates) {
    if (cov %in% missing) {
      cat(sprintf("MISSING : '%s'\n", cov))
    } else {
      cat(sprintf("OK      : %-30s  type=%-10s  NAs=%d\n",
                  cov, class(cd[[cov]]), sum(is.na(cd[[cov]]))))
    }
  }
  if (length(missing) > 0) {
    cat("\nAvailable colData columns:\n")
    cat(paste(" ", colnames(cd), collapse = "\n"), "\n\n")
    stop("Missing covariates: ", paste(missing, collapse = ", "),
         "\nFix the covariates vector before running the pipeline.")
  }
  cat("Done.\n\n")
}

# Fits one gene's marginal model in a forked worker (via mcparallel() in
# run_marginal_stage()) and caches the result to `outdir`. Not exported;
# called only from run_marginal_stage().
#' @noRd
.run_one_gene <- function(g, idx, sce, outdir,
                          covariates, mu_formula, sigma_formula,
                          celltype, family_use, timeout_sec) {
  file_out <- file.path(outdir, paste0(idx, "_", g, ".rds"))
  if (file.exists(file_out)) {
    res <- tryCatch(readRDS(file_out), error = function(e) NULL)
    mi  <- if (!is.null(res)) res$model_info else list(family_resolved=NA, engine=NA, link=NA, theta=NA, nb_fallback_reason=NA)
    tl  <- if (!is.null(res)) res$time_log   else list(total=NA)
    mu  <- if (!is.null(res) && !is.null(res$mu_hat)) mean(res$mu_hat, na.rm=TRUE) else NA_real_
    return(list(ok=TRUE, skipped=TRUE,
                family  = mi$family_resolved,
                engine  = mi$engine,
                link    = mi$link,
                theta   = mi$theta,
                mu      = mu,
                fallback_reason = mi$nb_fallback_reason,
                elapsed = as.numeric(tl$total, units="secs")))
  }

  sce_g <- sce[g, , drop = FALSE]
  t0 <- Sys.time()
  result <- tryCatch({
    res <- R.utils::withTimeout(
      sc3workflow_no_copula(
        sce             = sce_g,
        celltype        = celltype,
        other_covariates = covariates,
        mu_formula      = mu_formula,
        sigma_formula   = sigma_formula,
        family_use      = family_use,
        n_cores         = 1L
      ),
      timeout   = timeout_sec,
      onTimeout = "error"
    )
    saveRDS(res, file_out)
    mi <- res$model_info
    tl <- res$time_log
    mu <- if (!is.null(res$mu_hat)) mean(res$mu_hat, na.rm=TRUE) else NA_real_
    rm(res); gc()
    list(ok=TRUE, skipped=FALSE,
         family  = mi$family_resolved,
         engine  = mi$engine,
         link    = mi$link,
         theta   = mi$theta,
         mu      = mu,
         fallback_reason = mi$nb_fallback_reason,
         elapsed = as.numeric(tl$total, units="secs"))
  }, error = function(e) {
    call_str <- tryCatch(deparse(conditionCall(e)), error = function(...) "")
    msg <- sprintf("[FAILED] %s: %s  (call: %s)", g, e$message, call_str)
    cat(paste0(msg, "\n"), file = file.path(outdir, "failed_genes.txt"), append = TRUE)
    elapsed <- as.numeric(Sys.time() - t0, units = "secs")
    list(ok=FALSE, skipped=FALSE,
         family=NA, engine=NA, link=NA, theta=NA, mu=NA, fallback_reason=NA, elapsed=elapsed,
         error_msg=e$message)
  })
  rm(sce_g); gc()

  cat(sprintf("%s\t%.1fs\t%s\n", g, result$elapsed,
              if (result$ok) "OK" else "FAILED"),
      file = file.path(outdir, "timing.txt"), append = TRUE)
  result
}

# RNG streams for the forked gene fits. Each gene is fit in its own
# forked worker via parallel::mcparallel()/mccollect() below, and the
# randomized PIT jitter inside sc3workflow_no_copula() (the runif() draw
# behind U_randomized/Z_randomized) uses the RNG state the child receives
# at fork time.
# - fix_rng_streams = TRUE (default) sets RNGkind("L'Ecuyer-CMRG") for the
#   duration of the call (restored via on.exit()), so every forked gene
#   fit gets its own independent substream. The parallel package keeps
#   the stream it hands to forked children across calls within one R
#   session, and set.seed() alone does not reset it; set.seed() followed
#   by parallel::mc.reset.stream() does. The `seed` argument does exactly
#   this, which makes the output bit-identical across repeated calls and
#   across n_cores (gene i is always the i-th mcparallel() call). In a
#   fresh R session, set.seed() before the call is also sufficient.
# - fix_rng_streams = FALSE leaves the RNG kind unchanged, as in earlier
#   versions. Under R's default Mersenne-Twister, mcparallel() re-seeds
#   each forked child from the time and process ID, so the jitter cannot
#   be reproduced by any seed.
# Without a seed, repeated runs differ slightly because of the jitter.
# The effect on the shrunk estimates is small; see the package vignette
# and the associated paper.

#' Fit Stage 1's marginal model for every gene
#'
#' Fits a covariate-adjusted marginal model (negative binomial by
#' default) to each gene independently, in its own forked worker, via
#' `sc3workflow_no_copula()` (internal), and converts each fitted model's
#' probability-integral transform into Gaussian-copula scores. Results
#' are cached per gene under `comb_path/gene_fits/`, then assembled into
#' cell x gene CSV matrices written to `comb_path`.
#'
#' @param sce A `SingleCellExperiment`, genes x cells, with all of
#'   `covariates` present in `colData(sce)`. Should already be subset
#'   to a single cell type: POPCORN's donor-level correlation
#'   structure is defined within one cell type, and `celltype` below
#'   does not perform that subsetting for you.
#' @param comb_path Root output directory (created if absent). Per-gene
#'   fits and output CSVs under here are cached by filename ("already
#'   exists, skip"); re-calling with a different `sce` (different cells
#'   or genes) at the same `comb_path` errors rather than silently
#'   reusing mismatched cached output -- use a fresh `comb_path`, or
#'   `unlink(comb_path, recursive = TRUE)` first.
#' @param donor_col `colData(sce)` column holding donor identity. No
#'   default -- POPCORN estimates donor-*specific* correlations, so
#'   which column identifies donors is a modeling decision every
#'   caller must make explicitly, not something to silently fall back
#'   on. Kept in the saved metadata even when not one of `covariates`,
#'   so that [run_donor_split()] can find it. We recommend also
#'   including it in `covariates`/`mu_formula` (see there), so the
#'   mean model itself adjusts for between-donor compositional
#'   differences, on top of Stage 2's separate between-donor
#'   heterogeneity estimation.
#' @param gene_names Character vector of genes to assemble into the
#'   output matrices; must be a subset of `rownames(sce)` (`sce` may
#'   hold more genes than this -- unused ones are simply ignored). If
#'   `NULL` (default), every gene in `rownames(sce)` is used.
#' @param covariates Character vector of `colData(sce)` columns to pass
#'   to the marginal model. Default `"log1p_total_counts"` alone; we
#'   recommend also adding `donor_col` (e.g.
#'   `c("log1p_total_counts", "Donor")`) and any other
#'   technical/compositional covariates your data needs, and listing
#'   them in `mu_formula` too.
#' @param mu_formula GAMLSS mu formula (right-hand side only), e.g.
#'   `"log1p_total_counts + Donor"`, or
#'   `"log1p_total_counts + percent_mito + Sex + Donor"`. Default
#'   `"log1p_total_counts"` alone, matching the default `covariates`.
#' @param sigma_formula GAMLSS sigma formula. Default `"1"` (no
#'   covariate effect on dispersion), used throughout the associated
#'   paper.
#' @param family_use `"nb"`, `"poisson"`, `"gaussian"`, or `"auto"`.
#' @param celltype `colData(sce)` column holding a cell-type-like
#'   label. `scDesign3::construct_data()` requires at least one of
#'   `celltype`/`pseudotime`/`spatial` to be non-`NULL` (this package
#'   always passes the latter two as `NULL`). Default `"cell_type"`;
#'   if `sce` has no column by this name, one is added automatically
#'   with a single constant level, so most callers can ignore this
#'   argument entirely. If you already have a real column here (e.g.
#'   from a dataset with a genuine cell-type label), it is used as
#'   given -- but note it still does not subset `sce` by cell type;
#'   do that yourself beforehand if `sce` mixes more than one.
#' @param n_cores Number of parallel workers for gene fitting.
#' @param timeout_sec Per-gene timeout, in seconds.
#' @param blas_threads_per_gene Caps BLAS/OpenMP threads per gene fit, so
#'   that a single `gam()`/`gamlss()` call cannot saturate every core on
#'   the machine regardless of `n_cores`.
#' @param var_names Which per-gene quantities to assemble into output
#'   CSV matrices.
#' @param fix_rng_streams If `TRUE` (default), sets
#'   `RNGkind("L'Ecuyer-CMRG")` for the duration of the call (restored on
#'   exit) so that each forked gene fit draws from its own independent RNG
#'   substream. Required for `seed`. `FALSE` leaves the RNG kind
#'   unchanged, as in earlier versions; the randomized-PIT jitter then
#'   cannot be reproduced by any seed.
#' @param seed Optional integer. If given, the call runs `set.seed(seed)`
#'   followed by [parallel::mc.reset.stream()] under L'Ecuyer-CMRG before
#'   fitting, so the randomized-PIT jitter, and hence every output, is
#'   bit-reproducible, including across repeated calls in one R session
#'   and across `n_cores`. Like `set.seed()`, this changes the global RNG
#'   state. Default `NULL` leaves the RNG untouched.
#'
#' @section Run-to-run variability:
#' POPCORN's marginal fit uses a randomized probability-integral
#' transform, so two unseeded runs on identical data give slightly
#' different raw correlations. This is a deliberate statistical design
#' choice, not a bug, and its magnitude is small relative to the
#' heterogeneous edges reported in the associated paper. Pass `seed` for
#' an exactly reproducible run. A bare `set.seed()` before the call is
#' enough only in a fresh R session, because the parallel package keeps
#' the random-number stream it hands to forked workers across calls and
#' `set.seed()` does not reset it.
#'
#' @return Invisibly, a list with `comb_path` and `gene_path`.
#' @export
run_marginal_stage <- function(
  sce,
  comb_path,
  donor_col,
  gene_names    = NULL,
  covariates    = c("log1p_total_counts"),
  mu_formula    = "log1p_total_counts",
  sigma_formula = "1",
  family_use    = "nb",
  celltype      = "cell_type",
  n_cores       = 20L,
  timeout_sec   = 3600L,
  blas_threads_per_gene = 2L,
  var_names     = c("U_randomized_raw", "U_randomized",
                    "Z_randomized_raw", "Z_randomized",
                    "mu_hat", "sigma_hat"),
  fix_rng_streams = TRUE,
  seed          = NULL
) {
  if (!is.null(seed) && !fix_rng_streams) {
    stop("`seed` requires fix_rng_streams = TRUE: without L'Ecuyer-CMRG ",
         "streams, forked gene fits are re-seeded from time and PID.",
         call. = FALSE)
  }
  if (is.null(gene_names)) gene_names <- rownames(sce)

  # scDesign3::construct_data() requires at least one of
  # celltype/pseudotime/spatial to be non-NULL; this package always
  # passes the latter two as NULL (single-cell-type use only), so
  # `celltype` is the only one ever supplied. Most callers have no
  # real cell-type column to give it (sce is already single-cell-type
  # by the time it gets here), so add one with a single constant level
  # if it isn't already there, rather than requiring every caller to
  # do this themselves.
  if (!celltype %in% colnames(SummarizedExperiment::colData(sce))) {
    SummarizedExperiment::colData(sce)[[celltype]] <- factor("all")
  }

  if (fix_rng_streams) {
    old_rngkind <- RNGkind()
    on.exit(RNGkind(old_rngkind[1]), add = TRUE)
    RNGkind("L'Ecuyer-CMRG")
    if (!is.null(seed)) {
      # set.seed() alone does not reset parallel's forked-worker stream
      # within an R session; mc.reset.stream() re-derives it from the
      # freshly seeded .Random.seed (see REPRODUCIBILITY note above).
      set.seed(seed)
      parallel::mc.reset.stream()
    }
  }

  # gam()'s internal matrix operations let OpenBLAS spin up threads across
  # nearly all CPU cores automatically, independent of n_cores (gene-level
  # parallelism): even with n_cores=1 (fitting one gene at a time,
  # sequentially), a single gam() call can still saturate every core on the
  # machine. This caps each gene's fit at blas_threads_per_gene threads, so
  # the whole machine isn't saturated regardless of how n_cores is set.
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(blas_threads_per_gene)
    RhpcBLASctl::omp_set_num_threads(blas_threads_per_gene)
  }
  Sys.setenv(
    OMP_NUM_THREADS      = as.character(blas_threads_per_gene),
    OPENBLAS_NUM_THREADS = as.character(blas_threads_per_gene),
    MKL_NUM_THREADS      = as.character(blas_threads_per_gene)
  )

  gene_path <- file.path(comb_path, "gene_fits")
  dir.create(comb_path,  showWarnings = FALSE, recursive = TRUE)
  dir.create(gene_path,  showWarnings = FALSE, recursive = TRUE)

  # Per-gene fits and the U/Z_randomized CSVs below are cached by filename
  # only ("already exists, skip"), with no built-in check that a pre-existing
  # comb_path was written for the same dataset. Re-running with a different
  # sce (different cells and/or genes) at the same comb_path would silently
  # mix stale cached rows with the new cell_covariates written below, which
  # corrupts Stage 1b/2's downstream indexing instead of erroring cleanly.
  # Catch that mismatch here, before any stale file gets reused.
  meta_path <- file.path(comb_path, "metadata.rds")
  if (file.exists(meta_path)) {
    old_meta <- readRDS(meta_path)
    if (!identical(old_meta$cells, colnames(sce)) || !identical(old_meta$genes, rownames(sce))) {
      stop(sprintf(paste(
        "run_marginal_stage: comb_path '%s' already has cached Stage 1 output",
        "for a different dataset (%d cells, %d genes), but the current sce",
        "has %d cells, %d genes. Reusing it would silently mix stale and fresh",
        "rows rather than erroring cleanly downstream.",
        "Use a different comb_path, or clear this one first:",
        "  unlink(\"%s\", recursive = TRUE)",
        sep = "\n"
      ), comb_path, length(old_meta$cells), length(old_meta$genes),
         ncol(sce), nrow(sce), comb_path))
    }
  }

  check_covariates(sce, covariates)

  # --- fit one gene per core ---
  genes <- rownames(sce)
  message(sprintf("[Stage 1: marginal fit] Fitting %d genes with %d cores (family=%s) ...",
                  length(genes), n_cores, family_use))

  # Guards against mcparallel/mccollect returning a malformed result when a
  # forked worker is killed/interrupted before it can serialize its return
  # value (e.g. IRkernel's ZMQ socket receiving EINTR around a fork boundary
  # inside Jupyter). In that case mccollect hands back a bare atomic vector
  # (often class "try-error") instead of our expected list -- treat that gene
  # as failed and move on rather than letting `r$ok` crash the whole run.
  .sanitize_result <- function(r, g, outdir) {
    if (is.list(r) && !is.null(r$ok)) return(r)
    reason <- if (inherits(r, "try-error")) {
      paste("worker process crashed/interrupted:", conditionMessage(attr(r, "condition")))
    } else {
      "worker process returned no result (crashed or interrupted before completion)"
    }
    cat(sprintf("[FAILED] %s: %s\n", g, reason),
        file = file.path(outdir, "failed_genes.txt"), append = TRUE)
    list(ok = FALSE, skipped = FALSE,
         family = NA, engine = NA, link = NA, theta = NA, mu = NA, fallback_reason = NA, elapsed = NA,
         error_msg = reason)
  }

  .print_gene_result <- function(i, n, g, r) {
    if (isTRUE(r$ok)) {
      tag <- if (isTRUE(r$skipped)) "SKIP" else "OK  "
      cat(sprintf("[%s] [%3d/%3d] %-22s  family=%-14s  engine=%-8s  link=%-10s  mu=%-10s  theta=%-10s  %.1fs\n",
                  tag, i, n, g,
                  if (is.na(r$family)) "?" else r$family,
                  if (is.na(r$engine)) "?" else r$engine,
                  if (is.na(r$link))   "?" else r$link,
                  if (is.null(r$mu)    || is.na(r$mu))    "?" else signif(r$mu, 4),
                  if (is.null(r$theta) || is.na(r$theta)) "?" else signif(r$theta, 4),
                  if (is.na(r$elapsed)) 0 else r$elapsed))
      if (!is.null(r$fallback_reason) && !is.na(r$fallback_reason)) {
        cat(sprintf("       -> NB fallback reason: %s\n", r$fallback_reason))
      }
    } else {
      cat(sprintf("[FAIL] [%3d/%3d] %-22s  %.1fs  ERROR: %s\n",
                  i, n, g,
                  if (is.na(r$elapsed)) 0 else r$elapsed,
                  if (!is.null(r$error_msg)) r$error_msg else "unknown"))
    }
    flush(stdout())
  }

  if (n_cores == 1L) {
    # mcparallel+mccollect: parent never calls gam(), so C-level mgcv workspace
    # never accumulates in the parent.  Each child forks from a clean parent,
    # fits one gene, saves the RDS, then exits -- freeing all leaked C memory.
    # Copy-on-write means the large SCE is not physically copied per child.
    for (i in seq_along(genes)) {
      job <- parallel::mcparallel(
        .run_one_gene(
          g             = genes[i],
          idx           = i,
          sce           = sce,
          outdir        = gene_path,
          covariates    = covariates,
          mu_formula    = mu_formula,
          sigma_formula = sigma_formula,
          celltype      = celltype,
          family_use    = family_use,
          timeout_sec   = timeout_sec
        )
      )
      r <- parallel::mccollect(job, wait = TRUE)[[1]]
      r <- .sanitize_result(r, genes[i], gene_path)
      .print_gene_result(i, length(genes), genes[i], r)
    }
  } else {
    # Worker-pool via mcparallel/mccollect(wait=FALSE): keeps up to n_cores jobs
    # in flight and prints each gene's result as soon as it finishes, instead
    # of blocking (silently) until the whole mclapply batch completes.
    n <- length(genes)
    pending <- list()          # job -> gene index
    next_i  <- 1L
    run_results <- vector("list", n)

    launch_next <- function() {
      i <- next_i
      next_i <<- next_i + 1L
      job <- parallel::mcparallel(
        .run_one_gene(
          g             = genes[i],
          idx           = i,
          sce           = sce,
          outdir        = gene_path,
          covariates    = covariates,
          mu_formula    = mu_formula,
          sigma_formula = sigma_formula,
          celltype      = celltype,
          family_use    = family_use,
          timeout_sec   = timeout_sec
        )
      )
      pending[[as.character(job$pid)]] <<- i
    }

    while (next_i <= n && length(pending) < n_cores) launch_next()

    while (length(pending) > 0) {
      done <- parallel::mccollect(wait = FALSE)
      for (pid_chr in names(done)) {
        if (!pid_chr %in% names(pending)) next
        i <- pending[[pid_chr]]
        r <- .sanitize_result(done[[pid_chr]], genes[i], gene_path)
        run_results[[i]] <- r
        .print_gene_result(i, n, genes[i], r)
        pending[[pid_chr]] <- NULL
        if (next_i <= n) launch_next()
      }
      if (length(pending) > 0) Sys.sleep(0.2)
    }
  }

  # --- record formulas and metadata ---
  writeLines(
    c(paste0("mu_formula: ",    mu_formula),
      paste0("sigma_formula: ", sigma_formula),
      paste0("family_use: ",    family_use)),
    file.path(comb_path, "covariates.txt")
  )

  # Always include donor_col in cell_covariates so run_donor_split can find it,
  # even when donor_col is not part of the model covariates.
  cd           <- SummarizedExperiment::colData(sce)
  save_cols    <- union(covariates, intersect(donor_col, colnames(cd)))
  saveRDS(
    list(
      cell_covariates = cd[, save_cols, drop = FALSE],
      mu_formula      = mu_formula,
      sigma_formula   = sigma_formula,
      family_use      = family_use,
      genes           = rownames(sce),
      cells           = colnames(sce)
    ),
    file.path(comb_path, "metadata.rds")
  )

  # --- assemble per-gene results into cell x gene matrices and export CSV ---
  meta     <- readRDS(file.path(comb_path, "metadata.rds"))
  cell_cov <- meta$cell_covariates
  ngene    <- length(gene_names)

  message(sprintf("[Stage 1: marginal fit] Assembling %d genes into CSV matrices ...", ngene))

  read_gene_fit <- function(i) {
    fn <- file.path(gene_path, paste0(i, "_", gene_names[i], ".rds"))
    if (!file.exists(fn)) return(NULL)
    res <- readRDS(fn)
    stats::setNames(lapply(var_names, function(v) unname(res[[v]])), var_names)
  }

  res_list <- if (n_cores == 1L) {
    lapply(seq_len(ngene), read_gene_fit)
  } else {
    parallel::mclapply(seq_len(ngene), read_gene_fit, mc.cores = n_cores)
  }

  failed_idx  <- which(vapply(res_list, is.null, logical(1)))
  if (length(failed_idx) > 0) {
    warning(sprintf("[Stage 1: marginal fit] %d gene(s) missing fit files and will be excluded: %s",
                    length(failed_idx),
                    paste(gene_names[failed_idx], collapse = ", ")))
    res_list   <- res_list[-failed_idx]
    gene_names <- gene_names[-failed_idx]
    ngene      <- length(gene_names)
  }

  ncell <- length(res_list[[1]][[var_names[1]]])

  for (v in var_names) {
    out_csv <- file.path(comb_path, paste0(v, ".csv"))
    if (file.exists(out_csv)) next
    M <- matrix(unlist(lapply(res_list, "[[", v)), nrow = ncell, ncol = ngene)
    colnames(M) <- gene_names
    rownames(M) <- rownames(cell_cov)
    utils::write.csv(M, out_csv, row.names = TRUE)
    message("  wrote ", basename(out_csv))
  }

  message("[Stage 1: marginal fit] Done.")
  invisible(list(comb_path = comb_path, gene_path = gene_path))
}

# Fits one gene's covariate-adjusted marginal model (GAMLSS, by default
# negative binomial) via scDesign3's construct_data()/fit_marginal()/
# extract_para(), falling back to Poisson then Gaussian if the requested
# family fails to fit, and converts the fitted model into both the
# standard and randomized probability-integral-transform (PIT) scores.
# Not exported; called only from .run_one_gene(), one gene at a time.
#' @noRd
sc3workflow_no_copula <- function(
    sce,
    celltype = "Cell_type",
    other_covariates = c("log_total_count", "percent.mt", "batch"),
    family_use = "nb",
    mu_formula = "log_total_count + percent.mt + batch",
    sigma_formula = "1",
    n_cores = 1,
    verbose = FALSE
) {
    stopifnot(inherits(sce, "SingleCellExperiment"))
    stopifnot(nrow(sce) == 1)  # single-gene call

    time_log <- list()

    ## Step 0 (pre-process): merge all-zero donor levels into reference donor
    ## prevents rank-deficient (linearly dependent) model matrix columns from
    ## causing the GAM fit to fail
    gene_g <- rownames(sce)[1]
    x_g    <- as.numeric(SingleCellExperiment::counts(sce)[gene_g, ])
    for (cov_nm in other_covariates) {
        cd_now  <- SingleCellExperiment::colData(sce)
        if (!cov_nm %in% colnames(cd_now)) next
        grp_vec <- cd_now[[cov_nm]]
        if (!is.factor(grp_vec)) next
        grp_sums  <- tapply(x_g, grp_vec, sum)
        zero_lvls <- names(grp_sums)[grp_sums == 0]
        if (length(zero_lvls) == 0 || length(zero_lvls) == length(grp_sums)) next
        ## use the non-zero donor with the highest counts as the reference
        ## level, and merge every all-zero donor into it
        ref_lvl   <- names(which.max(grp_sums))
        new_grp   <- as.character(grp_vec)
        new_grp[new_grp %in% zero_lvls] <- ref_lvl
        keep_lvls <- setdiff(levels(grp_vec), zero_lvls)
        SummarizedExperiment::colData(sce)[[cov_nm]] <- factor(new_grp, levels = keep_lvls)
        message(sprintf(
            "  [pre-fix] gene '%s': merged %d all-zero %s levels into '%s' (rank fix)",
            gene_g, length(zero_lvls), cov_nm, ref_lvl))
    }

    ## Step 1: construct_data ----
    if (verbose) message("Step 1: construct_data")
    t1 <- Sys.time()
    data_obj <- scDesign3::construct_data(
        sce = sce,
        assay_use = "counts",
        celltype = celltype,
        pseudotime = NULL,
        spatial = NULL,
        other_covariates = other_covariates,
        corr_by = "1"
    )
    time_log$construct_data <- Sys.time() - t1

    ## Step 2: fit_marginal ----
    if (verbose) message("Step 2: fit_marginal")
    t2 <- Sys.time()
    family_resolved <- family_use
    if (family_use == "auto") {
        ## try nb -> poisson -> gaussian in order and use the first one that
        ## succeeds (the downstream PIT step only supports these three families)
        families_to_try <- c("nb", "poisson", "gaussian")
        marginal <- NULL
        for (fam in families_to_try) {
            marginal <- tryCatch(
                scDesign3::fit_marginal(
                    data = data_obj,
                    family_use = fam,
                    mu_formula = mu_formula,
                    sigma_formula = sigma_formula,
                    n_cores = n_cores
                ),
                error = function(e) NULL
            )
            if (!is.null(marginal)) {
                family_resolved <- fam
                if (verbose) message("  auto-selected family: ", fam)
                break
            }
        }
        if (is.null(marginal))
            stop("family_use = 'auto': all of nb, poisson, gaussian failed for this gene")
    } else {
        marginal <- scDesign3::fit_marginal(
            data = data_obj,
            family_use = family_use,
            mu_formula = mu_formula,
            sigma_formula = sigma_formula,
            n_cores = n_cores
        )
    }
    time_log$fit_marginal <- Sys.time() - t2

    ## check whether fit_marginal succeeded (scDesign3 catches errors
    ## internally and returns fit=NULL instead of raising)
    ## on failure, automatically fall back to poisson / gaussian
    gene_check <- rownames(sce)[1]
    .check_fit_null <- function(mg, g) {
        fc <- mg[[g]]
        is.null(fc$fit) || identical(fc$fit, NA)
    }
    .print_sc3_log <- function(mg, g) {
        fc <- mg[[g]]
        msgs <- character(0)
        if (!is.null(fc$warning) && length(fc$warning) > 0) {
            for (log_entry in fc$warning) {
                msg <- sprintf("[scDesign3 internal] %s [%s]: %s",
                                log_entry$function_name,
                                log_entry$type,
                                log_entry$message)
                message(msg)
                msgs <- c(msgs, msg)
            }
        }
        msgs
    }
    nb_fallback_reason <- NA_character_
    if (.check_fit_null(marginal, gene_check)) {
        nb_fallback_reason <- paste(.print_sc3_log(marginal, gene_check), collapse = " | ")
        if (nb_fallback_reason == "") nb_fallback_reason <- "NB fit returned NULL (no internal warning captured)"
        ## fallback attempt: poisson / gaussian
        fallback_done <- FALSE
        for (fb_fam in c("poisson", "gaussian")) {
            marginal_fb <- tryCatch(
                scDesign3::fit_marginal(
                    data = data_obj, family_use = fb_fam,
                    mu_formula = mu_formula,
                    sigma_formula = sigma_formula,
                    n_cores = n_cores
                ),
                error = function(e) NULL
            )
            if (!is.null(marginal_fb) && !.check_fit_null(marginal_fb, gene_check)) {
                marginal        <- marginal_fb
                family_resolved <- fb_fam
                message(sprintf("[fallback] gene '%s': NB failed -> using %s",
                                gene_check, fb_fam))
                fallback_done <- TRUE
                break
            }
        }
        if (!fallback_done) {
            stop(sprintf("fit_marginal failed for gene '%s': NB, poisson and gaussian all returned NULL fit",
                         gene_check))
        }
    }

    ## Step 3: extract_para ----
    if (verbose) message("Step 3: extract_para")
    t3 <- Sys.time()
    para <- scDesign3::extract_para(
        sce = sce,
        marginal_list = marginal,
        n_cores = n_cores,
        new_covariate = data_obj$new_covariate,
        data = data_obj$dat,
        family_use = family_resolved
    )
    time_log$extract_para <- Sys.time() - t3
    time_log$total <- sum(unlist(time_log))

    ## -----------------------------
    ## Step 4: identify fitted model
    ## -----------------------------
    gene <- rownames(sce)[1]
    fit_obj <- marginal[[gene]]$fit
    engine_used <- class(fit_obj)[1]   # "gamlss", "gam", etc.

    ## default NA
    fitted_family <- NA_character_
    fitted_link   <- NA_character_

    if (engine_used == "gam") {
        ## mgcv::gam
        fitted_family <- fit_obj$family$family
        fitted_link   <- fit_obj$family$link
    }

    if (engine_used == "gamlss") {
        ## gamlss
        fitted_family <- fit_obj$family[1]
        fitted_link   <- fit_obj$mu.link
    }

    ## Step 5: extract Y / mu / sigma ----
    y <- data_obj$count_mat[, 1]
    mu_hat <- para$mean_mat[, 1]
    sigma_hat <- if (!is.null(para$sigma_mat)) para$sigma_mat[, 1] else NULL

    ## -----------------------------------
    ## Step 6: compute U (standard <= , standard < , randomized PIT)
    ## -----------------------------------
    U_standard_raw     <- rep(NA_real_, length(y))  # P(X <= x)
    U_standard_lt_raw  <- rep(NA_real_, length(y))  # P(X <  x)
    U_randomized_raw   <- rep(NA_real_, length(y))
    eps <- 1e-6

    if (engine_used == "gam" && grepl("Negative Binomial", fitted_family)) {
        theta_msg <- fit_obj$family$getTheta(TRUE)
        message(
            "engine_used: ", engine_used,
            " | family: Negative Binomial",
            " | mu_mean: ", signif(mean(mu_hat), 4),
            " | mu_range: [", signif(min(mu_hat), 4), ", ", signif(max(mu_hat), 4), "]",
            " | theta: ", signif(theta_msg, 4)
        )
    } else {
        message(
            "engine_used: ", engine_used,
            " | family: ", fitted_family,
            " | mu_mean: ", signif(mean(mu_hat), 4)
        )
    }
    if (engine_used == "gamlss" && fitted_family %in% c("NB", "NBI", "NBII")) {
        ## gamlss NB (discrete): gamlss.dist provides pNBI / pNBII, there is
        ## no plain pNB
        pNB_fun <- if (fitted_family == "NBII") gamlss.dist::pNBII else gamlss.dist::pNBI
        Fy  <- pNB_fun(q = y,     mu = mu_hat, sigma = sigma_hat)     # P(X <= x)
        Fy1 <- pNB_fun(q = y - 1, mu = mu_hat, sigma = sigma_hat)     # P(X <  x)
        Fy2 <- pNB_fun(q = y - 2, mu = mu_hat, sigma = sigma_hat)     # P(X <= y-2), used for the lt randomized interval (Fy2, Fy1]

        U_standard_raw     <- Fy
        U_standard_lt_raw  <- Fy1
        U_randomized_raw  <- Fy1 + runif(length(y)) * (Fy - Fy1)
    } else if (engine_used == "gam" && fitted_family == "poisson") {
        ## Poisson (discrete)
        Fy  <- ppois(y,     lambda = mu_hat)
        Fy1 <- ppois(y - 1, lambda = mu_hat)
        Fy2 <- ppois(y - 2, lambda = mu_hat)

        U_standard_raw     <- Fy
        U_standard_lt_raw  <- Fy1
        U_randomized_raw   <- Fy1 + runif(length(y)) * (Fy - Fy1)
    } else if (engine_used == "gam" && fitted_family == "gaussian") {
        ## Gaussian (continuous)
        sd_hat <- if (!is.null(sigma_hat)) sigma_hat else sd(y, na.rm = TRUE)

        Fy <- pnorm(y, mean = mu_hat, sd = sd_hat)
        U_standard_raw     <- Fy
        U_standard_lt_raw  <- Fy      # continuous distribution: P(X < x) = P(X <= x)
        U_randomized_raw   <- Fy     # randomized PIT not needed

    } else if (engine_used == "gam" && grepl("Negative Binomial", fitted_family)) {
        ## mgcv NB (discrete)
        theta <- fit_obj$family$getTheta(TRUE)

        Fy  <- pnbinom(y,     mu = mu_hat, size = theta)
        Fy1 <- pnbinom(y - 1, mu = mu_hat, size = theta)
        Fy2 <- pnbinom(y - 2, mu = mu_hat, size = theta)

        U_standard_raw     <- Fy
        U_standard_lt_raw  <- Fy1
        U_randomized_raw   <- Fy1 + runif(length(y)) * (Fy - Fy1)
    } else {
        stop(
            sprintf(
                "Unsupported PIT combination: engine=%s, family=%s",
                engine_used, fitted_family
            )
        )
    }

    ## ----------------------------
    ## numerical safeguard (applied consistently to all four U variants)
    ## ----------------------------
    Z_standard_raw     <- qnorm(U_standard_raw)
    Z_standard_lt_raw  <- qnorm(U_standard_lt_raw)
    Z_randomized_raw   <- qnorm(U_randomized_raw)

    U_standard     <- pmin(pmax(U_standard_raw,     eps), 1 - eps)
    U_standard_lt  <- pmin(pmax(U_standard_lt_raw,  eps), 1 - eps)
    U_randomized   <- pmin(pmax(U_randomized_raw,   eps), 1 - eps)

    Z_standard     <- qnorm(U_standard)
    Z_standard_lt  <- qnorm(U_standard_lt)
    Z_randomized   <- qnorm(U_randomized)

    ## -----------------------------
    ## Step 7: return
    ## -----------------------------
    return(list(
        ## raw data
        Y = y,

        ## PIT (two variants)
        U_standard_raw = U_standard_raw,
        U_standard_lt_raw = U_standard_lt_raw,
        U_randomized_raw = U_randomized_raw,

        U_standard   = U_standard,
        U_standard_lt = U_standard_lt,
        U_randomized = U_randomized,

        Z_standard_raw   = Z_standard_raw,
        Z_standard_lt_raw = Z_standard_lt_raw,
        Z_randomized_raw = Z_randomized_raw,

        Z_standard   = Z_standard,
        Z_standard_lt = Z_standard_lt,
        Z_randomized = Z_randomized,

        ## parameters
        mu_hat    = mu_hat,
        sigma_hat = sigma_hat,

        ## model info (key)
        model_info = list(
            engine = engine_used,
            family_used = fitted_family,
            family_init = family_use,
            family_resolved = family_resolved,  # the family actually selected when family_use = "auto"
            link   = fitted_link,
            theta  = if (exists("theta")) theta else NA_real_,
            nb_fallback_reason = if (exists("nb_fallback_reason")) nb_fallback_reason else NA_character_,
            pit    = c("standard", "randomized") #, "standard_lt")
        ),

        ## timing
        time_log = time_log
    ))
}
