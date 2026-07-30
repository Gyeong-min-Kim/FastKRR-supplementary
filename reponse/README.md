# Response-letter code for *FastKRR*

This directory contains the simulation code and precomputed results used to
produce **Table S1** in our response letter to the editor and reviewers for
**"FastKRR: An R Package for Efficient Kernel Ridge Regression with
RcppArmadillo"** (Article ID: 2025-144).

Table S1 was added specifically to address a comment from Reviewer
(1-review-5.txt), who asked us to compare the REML-based tuning used by the
`KRMM` package with `FastKRR`'s tuning options (both the CVST-based
`"fastCV"` procedure and `FastKRR`'s own `"REML"` option), using both the
exact method and the low-rank approximations. This comparison is reported
only in the response letter and is not part of the numbered tables in the
revised manuscript.

## Files

- `tables1.R`

  Self-contained script that reproduces Table S1. It:
  1. Generates a single training set (\(n = 1{,}000\), \(d = 3\)) and a held-out
     evaluation set (\(n_z = 1{,}200\)) from the same data-generating process used
     in Section 6 of the manuscript.
  2. Benchmarks computation time (`microbenchmark`, 100 replications) for
     `KRMM::Kernel_Ridge_MM` and for `fastkrr()` under five configurations:
     exact method with `"REML"`, exact method with `"fastCV"`, and the
     Nyström, pivoted Cholesky, and RFF approximations with `"fastCV"`.
  3. Repeats data generation, model fitting, and out-of-sample prediction
     over 100 replications to estimate predictive accuracy (MSE and MAE, with
     standard errors) for the same five configurations plus `KRMM`.
  4. Combines the timing and prediction summaries into a single formatted
     table (`make_tbl()`, printed via `xtable`), corresponding to Table S1 in
     the response letter.

- `output/reml_vs_cvst_exact_time.Rdata`

  Precomputed `microbenchmark` timing results (object `exact_time`) used to
  produce the "Time" columns of Table S1.

- `output/reml_vs_cvst_exact_pred.Rdata`

  Precomputed prediction-error results across the 100 replications (object
  `compare_exact_pred`) used to produce the "Prediction" (MSE/MAE) columns of
  Table S1.

## Running the script

1. Install the required packages:

```r
install.packages(c("dplyr", "xtable", "ggplot2", "microbenchmark", "KRMM"))
```

2. Install `FastKRR` (CRAN v0.1.2 or later):

```r
install.packages("FastKRR")
```

3. Run the script from this directory:

```r
source("tables1.R")
```

Running the full script (100 timing replications plus 100 prediction
replications) takes a few minutes. It will (re)create the `output/` folder and
overwrite `output/reml_vs_cvst_exact_time.Rdata` and
`output/reml_vs_cvst_exact_pred.Rdata` with newly generated results, then print
the final table (`tbl_exact`).

## Reproducibility notes

The precomputed `.Rdata` files in `output/` reproduce the values reported as Table S1 in
the response letter. As with the other timing results in this repository,
exact runtimes depend on the local hardware and system load and are not
expected to match exactly when the script is rerun, although the relative
ordering of methods (`KRMM` far slower than any `FastKRR` option; `"REML"`
faster than `"fastCV"` for the exact method) should remain consistent.
