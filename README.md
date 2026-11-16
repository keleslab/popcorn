# popcorn

Population-Scale Personalized Coexpression Network Estimation.

`popcorn` estimates reliable, donor-specific gene co-expression
correlations from population-scale single-cell RNA-seq data. It works
in two stages: Stage 1 fits a covariate-adjusted marginal model to
correct technical and compositional noise within each donor's own
cells; Stage 2 applies a random-effects variance decomposition that
separates genuine between-donor heterogeneity from sampling noise at
each gene pair, adaptively shrinking noisy donor-level correlation
estimates toward the global correlation using information borrowed
jointly across donors.

## Installation

popcorn requires R (>= 4.3) on Linux or macOS (Stage 1 fits genes in
forked worker processes, which are not available on Windows). Install
from GitHub; Bioconductor dependencies (`SingleCellExperiment`,
`SummarizedExperiment`, `qvalue`) and `scDesign3` are installed along
with it:

```r
# install.packages(c("remotes", "BiocManager"))
options(repos = BiocManager::repositories())
remotes::install_github("keleslab/popcorn")
```

## Usage

```r
library(popcorn)

result <- run_coexp_pipeline(
  sce           = sce,            # a SingleCellExperiment, genes x cells
  comb_path     = "path/to/output",
  gene_names    = rownames(sce),
  covariates    = c("log1p_total_counts", "percent_mito", "Donor"),
  mu_formula    = "log1p_total_counts + percent_mito + Donor",
  sigma_formula = "1",
  donor_col     = "Donor"
)

head(result$edge_df)
```

See `vignette("popcorn-intro", package = "popcorn")` for a full worked
example on a simulated toy dataset, including how to interpret the
output (`tau2_e`, `p_Q_e`, the raw vs. shrunk per-donor correlations).

## Reproducibility

Stage 1 uses a randomized probability-integral transform. Pass `seed`
to `run_coexp_pipeline()` for a bit-reproducible run. Code reproducing
the analyses in the associated paper is at
<https://github.com/keleslab/popcorn-paper> and uses popcorn 0.1.1.

## Authors

Shuchen Yan, Kris Sankaran, Sündüz Keleş.

## License

MIT; see [LICENSE](LICENSE).
