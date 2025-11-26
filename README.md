# Supplementary materials for *FastKRR*

This repository provides the supplementary materials for the paper
**“FastKRR: An R Package for Efficient Kernel Ridge Regression with RcppArmadillo”**,
including the full replication scripts and the precomputed `.RData` files used to
generate all tables and figures in the manuscript.


## Why precomputed RData files are provided

This repository contains the computational resources used in two major components
of the analysis presented in the paper.
The first is the **numerical study**, which involves large-scale
simulation experiments and comparisons of computation time between the exact and
approximate FastKRR estimators. Because the timing component naturally depends
on the computing environment (e.g., CPU architecture, parallelism, and OS
scheduling), the exact runtime values cannot be reproduced identically across
systems, even with fixed random seeds.
To ensure that the results reported in the manuscript can be recreated
immediately, we provide the precomputed `.RData` files used in the tables and figures.

Following the numerical study, the paper includes a **tidymodels integration
example** to illustrate how FastKRR can be used within a modern modeling
workflow. This example does not require any precomputed data; it can be run
directly from the provided scripts. Although the computations may take some
time, the results are fully reproducible.

To support both immediate reproducibility and full transparency, we provide:

- **full replication scripts**, covering both the simulation study and the
  tidymodels workflow, and  
- **precomputed `.RData` files** *only for the simulation study*, generated on
  the system used during manuscript preparation.

Users may run the entire pipeline on their own machines, or rely on the
precomputed files to reproduce the figures and tables without rerunning the full
set of simulations.


## Main contents

- `scripts/full_replication_code.R`  
  Full pipeline to reproduce:
  - all simulation results (runtime: ~3–4 hours),
  - the tidymodels workflow example,
  - benchmarking and comparison results.  
  Running this script will generate new `.RData` files in  
  **`scripts/output/`**.

- `supplement_materials/*.RData`  
  Precomputed results used directly in the paper for tables and figures.

- Additional helper scripts and documentation.



## Reproducibility notes

The precomputed `.RData` files reproduce exactly the numerical values used in the
paper.  
Runtime values may differ across systems, but the **relative performance,
comparative patterns, and substantive conclusions remain unchanged**.