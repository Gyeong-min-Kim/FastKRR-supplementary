# Code

The `scripts - first-revision/` folder includes `full_replication_code-first-rivision.R`, which provides 
the complete replication pipeline for all analyses reported in **Sections 5–7** 
of the manuscript.


# Running the replication scripts

1. Install the necessary packages:

```r
install.packages(c("dplyr", "tidyr", "ggplot2", "tidymodels",
                   "microbenchmark", "knitr, "kableExtra", "mgcv"
                   "xtable", "KRLS", "gKRLS", "kernlab", "FastKRR"))
```

2. The package FastKRR **v0.2.1** is available on CRAN: https://cran.r-project.org/package=FastKRR
  `install.packages("FastKRR")`


3. Load the package: `library(FastKRR)`


4. To reproduce the full results, run the replication scripts from the `scripts - first-revision/` directory


# Note
Running `full_replication_code-first-rivision.R` reproduces the results reported in the manuscript 
using the same analysis pipeline, and new `.RData` or `.Rdata` files will be created in 
the `scripts - first-revision/output/` folder.
The timing values in the replication `.RData` or `.Rdata` files may differ slightly from those
in the manuscript and can vary across systems due to hardware, OS scheduling,
or parallelization. However, the substantive results and overall patterns remain 
consistent, and the performance comparisons and conclusions are unaffected.
The exact `.RData` or `.Rdata` files used in the manuscript (not newly generated ones)
are available at: https://github.com/Gyeong-min-Kim/FastKRR-supplementary