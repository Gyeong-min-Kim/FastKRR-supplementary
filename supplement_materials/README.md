# Supplemental Materials

This folder contains `.RData` files that were generated on our reference
machine and are used to compile the manuscript exactly as submitted.
They ensure that all tables and figures in the paper reproduce
bit-for-bit without rerunning the full replication code.

- These files correspond to intermediate results from the numerical
  studies (e.g., `package_comp.RData`, `plot_n5000_d3_eps6.RData`).
- They are loaded by `FastKRR.pdf` during compilation.
- Their presence guarantees that the compiled manuscript matches the
  published version, even though some computations would otherwise take
  several hours.

Note: All of these `.RData` files can also be regenerated from scratch
by running `scripts/full_replication_code.R`. The regenerated files
will be saved in `scripts/output/` and may differ slightly in timing
values due to hardware and OS variation, but the substantive results and
conclusions remain the same.