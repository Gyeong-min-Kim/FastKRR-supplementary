# Supplementary materials for *FastKRR*

This repository provides the supplementary materials for the paper
**“FastKRR: An R Package for Efficient Kernel Ridge Regression with RcppArmadillo”**,
including replication scripts and precomputed `.Rdata` files used
to generate the tables and figures in the submission.

## Why precomputed RData files are provided

This repository contains the computational resources used in two major components
of the analysis presented in the paper.
The **numerical study** includes large-scale simulation experiments,
comparisons of computation time and predictive accuracy, evaluations of exact
and approximate KRR implementations, and comparisons with other R packages.
Because runtime measurements depend on the computing environment, including
CPU architecture, parallel configuration, operating system, and system
scheduling, exact computation times may not be reproduced identically across
machines, even when random seeds are fixed.
To allow the numerical values reported in the manuscript to be reproduced
without rerunning all time-consuming experiments, we provide the precomputed
`.Rdata` files used to generate the reported tables and figures.

The paper also includes a **tidymodels integration example** illustrating how
FastKRR can be incorporated into a modern modeling workflow. This example does
not require precomputed result files and can be run directly from the provided
replication scripts.


## Repository organization

The supplementary files are organized according to the stage at which each
analysis was introduced or updated.

- Files stored directly in `supplement_materials/` correspond to analyses
  included in the submission.

The manuscript may use files. Analyses retained
without modification continue to use the original result files.

## Main contents

- `scripts/full_replication_code.R`: Full pipeline to reproduce,
  - all simulation results (runtime: ~34–37 hours),
  - experiments with larger sample sizes and dimensions,
  - computation-time and predictive-accuracy evaluations,
  - the tidymodels workflow example,
  - benchmarking and comparison results.  
  - running this script will generate new `.Rdata` files in **`scripts/output/`**.

- `supplement_materials/*.Rdata`: Precomputed results used directly in the paper for tables and figures.

- `reponse/`: Simulation code and precomputed results for **Table S1** in the
  response letter to the editor and reviewers, comparing `KRMM`'s REML-based
  tuning with `FastKRR`'s `"fastCV"` and `"REML"` options.

- Additional helper scripts and documentation.

## Reproducibility notes

The precomputed result files reproduce the numerical values used in the
corresponding version of the manuscript.

Runtime values obtained by rerunning the experiments may vary across computing
environments. Consequently, exact timing values are not expected to be
identical across machines, although the main comparative patterns and
substantive conclusions should remain qualitatively similar.
