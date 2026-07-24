# First-revision result files for *FastKRR*

This directory contains the precomputed `.RData` files used for the numerical
experiments added or updated in the first revision of the paper
**“FastKRR: An R Package for Efficient Kernel Ridge Regression with RcppArmadillo.”**

The files were generated on the reference computing system used during the
revision and are provided so that the revised tables and figures can be
reproduced without rerunning the full set of computationally intensive
experiments.


## Why precomputed RData files are provided

The revised numerical study includes large-scale simulation experiments,
comparisons with external KRR-related packages, evaluations of computation time
and predictive accuracy, and comparisons of hyperparameter-selection
procedures.

Runtime measurements depend on the computing environment, including CPU
architecture, operating system, parallel configuration, system load, and
scheduling. Consequently, exact computation times may differ when the
experiments are rerun on another machine, even when the random seeds and model
settings are fixed.

The precomputed files in this directory contain the numerical results used to
generate the corresponding tables and figures in the first revised manuscript.


# File organization

### Revised Table 2: package-level computation-time comparison

- `table2_compare_time_packages_fastkrr(n1000d3).RData`

  Contains the package-level computation-time comparison among
  `FastKRR`, `KRLS`, `gKRLS`, and `kernlab` under fixed hyperparameter
  settings.

### Revised Table 3: FastKRR computation time and predictive accuracy

- `table3_compare_time_opt(all).Rdata`

  Contains computation-time results for the exact, Nyström, pivoted Cholesky,
  and random Fourier feature implementations across the revised simulation
  settings.

- `table3_compare_prediction (n5000d3).Rdata`
- `table3_compare_prediction (n5000d5).Rdata`
- `table3_compare_prediction (n5000d7).Rdata`
- `table3_compare_prediction (n10000d3).Rdata`
- `table3_compare_prediction (n10000d5).Rdata`
- `table3_compare_prediction (n10000d7).Rdata`

  Contain predictive-accuracy results for
  \(n \in \{5{,}000,10{,}000\}\) and \(d \in \{3,5,7\}\).

### Revised Table 4: approximation-method comparison

- `table4_compare_time_pkg(n10000d3_m111).Rdata`
- `table4_compare_time_pkg(n10000d3_m208).Rdata`

  Contain computation-time results for the package-level approximation
  comparison at approximation ranks \(m=111\) and \(m=208\).

- `table4_compare_prediction(n10000d3_m111).Rdata`
- `table4_compare_prediction(n10000d3_m208).Rdata`

  Contain predictive-accuracy results for the same approximation-rank
  settings.

  These analyses compare the sketched implementation in `gKRLS`, the Nyström
  implementation in `KRLS`, and the Nyström, pivoted Cholesky, and random
  Fourier feature approximations in `FastKRR`.

### Revised Table 5: hyperparameter-selection comparison

- `table5_compare_cv.Rdata`

  Contains computation-time results for ExactCV, FastCV, and REML within
  `FastKRR`.

- `table5_compare_cv_prediction (n2000d3).Rdata`

  Contains the corresponding predictive-accuracy results.

### Additional revised figure results

- `plot_n5000_d3_eps6.RData`

  Contains precomputed results used for a figure or graphical summary in the
  revised manuscript.


## Reproducibility notes

The files in this directory reproduce the numerical values reported in the
corresponding revised tables and figures.

Exact runtime values may vary across machines because of differences in
hardware, operating systems, parallel execution, and system load. Predictive
results may also exhibit minor floating-point or stochastic variation when
regenerated, but the main comparative patterns should remain qualitatively
similar.

Files corresponding to analyses retained unchanged from the original
submission remain in the parent `supplement_materials/` directory.