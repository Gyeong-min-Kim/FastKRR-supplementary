# Code

The `scripts/` folder includes `full_replication_code.R`, which provides 
the complete replication pipeline for all analyses reported in **Sections 5–7** 
of the manuscript.

Because the full simulations and benchmarking procedures are computationally intensive, 
running the script may require approximately **40 hours** (estimated from the recorded benchmark times; the actual runtime depends on hardware and parallelization settings).

# Running the replication script

1. Install the necessary packages (`FastKRR` **v0.2.2** is available on
   [CRAN](https://cran.r-project.org/package=FastKRR)):

```r
install.packages(c("dplyr", "tidyr", "ggplot2", "tidymodels", "tibble",
                   "microbenchmark", "knitr", "kableExtra", "mgcv",
                   "xtable", "KRLS", "gKRLS", "kernlab", "FastKRR"))
```

   The script calls `krls()` with `lambda_method`, `approx`, and `nystrom_m`, so a
   version of `KRLS` that supports these arguments is required (checked with
   `KRLS` 1.7-1).

   The `gKRLS` models are fitted using `mgcv::gam()`, with settings that differ by
   experiment (checked with `gKRLS` 1.0.4):

   - Table 2: `sketch_method = "none"`, `bandwidth = rho`, `standardize = "none"`,
     and fixed `sp = lambda`.
   - Table 4: `method = "REML"`, `sketch_method = "none"`, `bandwidth = rho`,
     `standardize = "none"`, and `remove_instability = FALSE`.
   - Table 5: `sketch_method = "subsampling"` with `sketch_size_raw = m`, using the
     default `gKRLS` bandwidth and standardization.

2. Run the script with `scripts/` as the working directory:

```r
setwd("scripts")  # from the repository root
source("full_replication_code.R")
```


# Note
Running `full_replication_code.R` reruns the analyses reported in the manuscript
using the same analysis pipeline, and new `.Rdata` files will be created in
the `scripts/output/` folder.
The timing values in the replication `.Rdata` files may differ from those
in the manuscript and can vary across systems due to hardware, OS scheduling,
or parallelization. However, the main comparative patterns and conclusions are
expected to remain qualitatively similar.
The exact `.Rdata` files used in the manuscript (not newly generated ones)
are available in `../supplement_materials/`.