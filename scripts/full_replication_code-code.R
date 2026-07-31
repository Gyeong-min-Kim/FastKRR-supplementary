rm(list = ls())
###### REQUIRED PACKAGES -------------------------------------------------------
library(FastKRR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidymodels)
library(microbenchmark)
library(knitr)
library(kableExtra)
library(KRLS)
library(gKRLS)
library(mgcv)
library(kernlab)
library(xtable)

start_time = Sys.time()
## 5. Structure and functionality of FastKRR package 
#generate data 
set.seed(1)
n = 1000 ; d = 1

X = matrix(runif(n*d, 0, 1), nrow = n, ncol = d)
y = as.vector(sin(2*pi*X^3) + rnorm(n, 0, 0.3))
data = data.frame(X = X, y = y)

X_new = matrix(seq(0, 1, length.out = 1200), ncol = 1)
y_new = as.vector(sin(2*pi*X_new^3) + rnorm(1200, 0, 0.3))
data_new = data.frame(X = X_new, y = y_new)


## 5.1.Main functions
#example package key arguments and options
fit_exact = fastkrr(data = data, response = "y", kernel = "gaussian", 
                    opt = "exact", selection_method = "exactCV")
fitted_exact = predict(fit_exact)
pred_exact = predict(fit_exact, newdata = X_new)

FastKRR::error(fit_exact)
FastKRR::error(fit_exact, data_new = data_new)

FastKRR::param(fit_exact)
summary(fit_exact)
plot(fit_exact)

#Comparison of kernel ridge regression fits using different computation methods
fit_exact   = fastkrr(data, response = "y", kernel = "gaussian", opt = "exact",  verbose = TRUE)
fit_nystrom = fastkrr(data, response = "y", kernel = "gaussian", opt = "nystrom", verbose = TRUE)
fit_pivoted = fastkrr(data, response = "y", kernel = "gaussian", opt = "pivoted", verbose = TRUE)
fit_rff     = fastkrr(data, response = "y", kernel = "gaussian", opt = "rff", verbose = TRUE)


pred_exact   = predict(fit_exact,   X_new)
pred_nystrom = predict(fit_nystrom, X_new)
pred_pivoted = predict(fit_pivoted, X_new)
pred_rff     = predict(fit_rff,     X_new)


df_preds = data.frame(
  X = X_new[,1],
  Exact   = pred_exact,
  Nystrom = pred_nystrom,
  Pivoted = pred_pivoted,
  RFF     = pred_rff) %>%
  pivot_longer(-X, names_to = "Method", values_to = "y_hat")

df_preds$Method = factor(
  df_preds$Method,
  levels = c("Exact", "Nystrom", "Pivoted", "RFF"),
  labels = c("(a) Exact", "(b) Nystrom", "(c) Pivoted Cholesky", "(d) RFF"))

df_train = data.frame(X = X[,1], y = y)

ggplot() +
  geom_point(data = df_train, aes(x = X, y = y),
             color = "grey70", alpha = 0.35, size = 0.6) +
  geom_line(data = df_preds, aes(x = X, y = y_hat, color = Method),
            linewidth = 0.6) +
  facet_wrap(~ Method, ncol = 2) +
  scale_color_manual(values = c(
    "(a) Exact"   = "#1f77b4",
    "(b) Nystrom" = "#d62728",
    "(c) Pivoted Cholesky" = "#2ca02c",
    "(d) RFF"     = "#9467bd"
  )) +
  labs(x = "X", y = "y") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")


fit_nystrom = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "nystrom")
fit_pivoted = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "pivoted")
fit_rff = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "rff")

## 5.2 Utility functions for kernel computations
K = make_kernel(X, kernel = "gaussian", rho = 2)
K_pred = make_kernel(X, X_new, kernel = "gaussian", rho = 2)


K = make_kernel(X, kernel = "laplace", rho = 2, n_threads = 4)
K_pred = make_kernel(X, X_new, kernel = "laplace", rho = 2, n_threads = 4)


## 5.3.Parallel computation support
K_nystrom = approx_kernel(X = X, rho = 1, opt = "nystrom", n_threads = 4)
K_pchol = approx_kernel(X = X, rho = 1, opt = "pivoted")
K_rff = approx_kernel(X = X, rho = 1, opt = "rff", n_threads = 4)

## 5.4.regularization parameter selection
#hyperparameter select examples
fit_default = fastkrr(data = data, response = "y", selection_method = "exactCV")

fit_fixed = fastkrr(data = data, response = "y", lambda = 1e-4)

fit_reml_range = fastkrr(data = data, response = "y", lambda = c(1e-8, 1e-2),
                         selection_method = "REML")

fit_grid = fastkrr(data = data, response = "y", lambda = 10^seq(-8, -2, length.out = 10),
                   selection_method = "fastCV")

## 6. Numerical study
#plot-funcion-------------------------------------------------------------------
time_plot = function(df, rm = 0, size = 15){
  df = as.data.frame(df)
  idx = df$expr %in% rm; df = df[!idx, ]
  stats = df %>% dplyr::group_by(expr) %>%
    dplyr::summarize(Q1 = stats::quantile(time/1e6, 0.25),
                     Q3 = stats::quantile(time/1e6, 0.75),
                     IQR = Q3 - Q1, .groups="drop")
  ymin = min(stats$Q1 - 1.5*stats$IQR); ymax = max(stats$Q3 + 1.5*stats$IQR)
  ggplot2::ggplot(df, ggplot2::aes(x = expr, y = time/1e6)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::coord_cartesian(ylim = c(ymin, ymax)) +
    ggplot2::labs(x = NULL, y = "Time (milliseconds)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = size))
}

#summary-functions--------------------------------------------------------------
summarize_package_comparison = function(
    benchmark_result,
    prediction_result,
    m_value = NULL,
    n_value = NULL,
    d_value = NULL
) {
  time_df = benchmark_result %>%
    as.data.frame() %>%
    mutate(expr = as.character(expr))
  
  if (!is.null(n_value) && !is.null(d_value)) {
    suffix = paste0("_n", n_value, "d", d_value, "$")
    time_df = time_df %>%
      filter(grepl(suffix, expr)) %>%
      mutate(
        Method = sub(suffix, "", expr),
        time_sec = time / 1e9)
  } else {
    time_df = time_df %>%
      mutate(
        Method = expr,
        time_sec = time / 1e9)
  }
  
  time_summary = time_df %>%
    group_by(Method) %>%
    summarise(
      Time_mean = mean(time_sec, na.rm = TRUE),
      
      Time_se = sd(time_sec, na.rm = TRUE) /
        sqrt(sum(!is.na(time_sec))),
      
      Time_median = median(time_sec, na.rm = TRUE),
      
      .groups = "drop"
    )
  
  prediction_long = purrr::imap_dfr(
    prediction_result,
    function(result_list, method_name) {
      
      purrr::map_dfr(
        result_list,
        function(result) {
          tibble(
            Method = method_name,
            MSE = mean(result^2, na.rm = TRUE),
            MAE = mean(abs(result), na.rm = TRUE))
        })})
  
  prediction_summary = prediction_long %>%
    group_by(Method) %>%
    summarise(
      MSE_mean = mean(MSE, na.rm = TRUE),
      
      MSE_se = sd(MSE, na.rm = TRUE) /
        sqrt(sum(!is.na(MSE))),
      
      MAE_mean = mean(MAE, na.rm = TRUE),
      
      MAE_se = sd(MAE, na.rm = TRUE) /
        sqrt(sum(!is.na(MAE))),
      
      .groups = "drop")
  
  result = time_summary %>%
    left_join(
      prediction_summary,
      by = "Method"
    ) %>%
    mutate(
      `Mean (SE)` = sprintf(
        "%.4f (%.4f)",
        Time_mean,
        Time_se
      ),
      
      Median = sprintf(
        "%.4f",
        Time_median
      ),
      
      `MSE (SE)` = sprintf(
        "%.4f (%.4f)",
        MSE_mean,
        MSE_se
      ),
      
      `MAE (SE)` = sprintf(
        "%.4f (%.4f)",
        MAE_mean,
        MAE_se
      )
    )
  
  
  if (!is.null(m_value)) {result$m = m_value}
  
  if (!is.null(n_value)) {result$n = n_value}
  
  if (!is.null(d_value)) {result$d = d_value}
  
  id_columns = c("Method", "m", "n", "d")
  
  id_columns = id_columns[id_columns %in% names(result)]
  
  result %>%
    select(
      all_of(id_columns),
      `Mean (SE)`,
      Median,
      `MSE (SE)`,
      `MAE (SE)`
    )
}

calc_metric = function(pkg_error_list) {
  mse = sapply(pkg_error_list, function(err) {
    err = as.vector(err)
    mean(err^2, na.rm = TRUE)
  })
  mae = sapply(pkg_error_list, function(err) {
    err = as.vector(err)
    mean(abs(err), na.rm = TRUE)
  })
  data.frame(
    MSE_mean = mean(mse, na.rm = TRUE),
    MSE_se   = sd(mse, na.rm = TRUE),
    MAE_mean = mean(mae, na.rm = TRUE),
    MAE_se   = sd(mae, na.rm = TRUE)
  )
}


## 6.1 Computation time with external packages ---------------------------------
#Computation-time-comparison ---------------------------------------------------
set.seed(1)
n = 1000
d = 3
X = matrix(runif(n*d,0, 1),n,d)
data = data.frame(X, y)

X1 = X[, 1]
X2 = X[, 2]
X3 = X[, 3]

y = as.vector(sin(2*pi*rowMeans(X)^3) + rnorm(n, 0, 0.1))

# hyper parameter
lambda  = 1e-4
rho = 1


# kernlab code

krr = function(X, y, lambda, rho){
  n_row = nrow(X)
  
  K = kernelMatrix(rbfdot(sigma=0.5), X)
  
  alpha_hat = solve(K + lambda * diag(n_row)) %*% y
  
  pred_y = t(alpha_hat) %*% K
  
  return(list("pred_y" = pred_y,
              "coefficient" = alpha_hat))
}

# time
krr_compare = microbenchmark::microbenchmark(
  KRLS = krls(X= X, y= as.vector(y), whichkernel = "gaussian",
              sigma = rho, derivative = FALSE, lambda = lambda,
              binary = FALSE, vcov = FALSE, approx = "none"),
  
  gKRLS = gam(y ~ s(X1, X2, X3, bs = "gKRLS", 
                    xt = gKRLS(
                      sketch_method = "none", 
                      bandwidth = rho,
                      standardize = "none"
                    )),
              sp = lambda),
  
  kernlab = krr(X, y, lambda = lambda, rho = rho),
  
  FastKRR = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "exact", 
                    rho = rho, lambda = lambda, verbose = FALSE),
  times = 100,
  unit = "ms"
)


save(krr_compare, file = "output/table2_compare_time_packages_fastkrr(n1000d3).Rdata")
df = as.data.frame(krr_compare)


summary_tbl = df %>%
  group_by(expr) %>%
  summarize(
    Mean   = mean(time)   / 1e9,
    SE     = sd(time) / sqrt(100) / 1e9,      
    Median = median(time) / 1e9,
    .groups = "drop"
  )


xt = xtable(summary_tbl, digits = 3,align = c("c", rep("c", ncol(summary_tbl))))
print(xt, include.rownames = FALSE, floating = FALSE)

## 6.2 Computation time and predictive accuracy in FastKRR
##hyper parameter
rho = 1

##select hyper parameter lambda each n and d
###n = 5000, d = 3
n = 5000
d = 3
X = matrix(runif(n*d,0, 1),n,d)
y = as.vector(sin(2*pi*rowMeans(X)^3) + rnorm(n, 0, 0.1))

data_n5000d3 = data.frame(X, y)

FastKRR_exact = fastkrr(data = data_n5000d3, response = "y", kernel = "gaussian", opt = "exact", 
                        rho = rho, selection_method = "REML", verbose = TRUE)
lambda_exact_n5000d3 = FastKRR_exact$lambda


FastKRR_nystrom = fastkrr(data = data_n5000d3, response = "y", kernel = "gaussian", opt = "nystrom", 
                          rho = rho, selection_method = "REML", verbose = TRUE)
lambda_nystrom_n5000d3 = FastKRR_nystrom$lambda


FastKRR_pivoted = fastkrr(data = data_n5000d3, response = "y", kernel = "gaussian", opt = "pivoted", 
                          rho = rho, selection_method = "REML", verbose = TRUE)
lambda_pivoted_n5000d3 = FastKRR_pivoted$lambda


FastKRR_rff = fastkrr(data = data_n5000d3, response = "y", kernel = "gaussian", opt = "rff", 
                      rho = rho, selection_method = "REML", verbose = TRUE)
lambda_rff_n5000d3 = FastKRR_rff$lambda


### n = 5000, d = 5
n = 5000
d = 5
X = matrix(runif(n*d,0, 1),n,d)
y = as.vector(sin(2*pi*rowMeans(X)^3) + rnorm(n, 0, 0.1))

data_n5000d5 = data.frame(X, y)

FastKRR_exact = fastkrr(data = data_n5000d5, response = "y", kernel = "gaussian", opt = "exact", 
                        rho = rho, selection_method = "REML", verbose = TRUE)
lambda_exact_n5000d5 = FastKRR_exact$lambda


FastKRR_nystrom = fastkrr(data = data_n5000d5, response = "y", kernel = "gaussian", opt = "nystrom", 
                          rho = rho, selection_method = "REML", verbose = TRUE)
lambda_nystrom_n5000d5 = FastKRR_nystrom$lambda


FastKRR_pivoted = fastkrr(data = data_n5000d5, response = "y", kernel = "gaussian", opt = "pivoted", 
                          rho = rho, selection_method = "REML", verbose = TRUE)
lambda_pivoted_n5000d5 = FastKRR_pivoted$lambda


FastKRR_rff = fastkrr(data = data_n5000d5, response = "y", kernel = "gaussian", opt = "rff", 
                      rho = rho, selection_method = "REML", verbose = TRUE)
lambda_rff_n5000d5 = FastKRR_rff$lambda


### n = 5000, d = 7
n = 5000
d = 7
X = matrix(runif(n*d,0, 1),n,d)
y = as.vector(sin(2*pi*rowMeans(X)^3) + rnorm(n, 0, 0.1))

data_n5000d7 = data.frame(X, y)

FastKRR_exact = fastkrr(data = data_n5000d7, response = "y", kernel = "gaussian", opt = "exact", 
                        rho = rho, selection_method = "REML", verbose = TRUE)
lambda_exact_n5000d7 = FastKRR_exact$lambda


FastKRR_nystrom = fastkrr(data = data_n5000d7, response = "y", kernel = "gaussian", opt = "nystrom", 
                          rho = rho, selection_method = "REML", verbose = TRUE)
lambda_nystrom_n5000d7 = FastKRR_nystrom$lambda


FastKRR_pivoted = fastkrr(data = data_n5000d7, response = "y", kernel = "gaussian", opt = "pivoted", 
                          rho = rho, selection_method = "REML", verbose = TRUE)
lambda_pivoted_n5000d7 = FastKRR_pivoted$lambda


FastKRR_rff = fastkrr(data = data_n5000d7, response = "y", kernel = "gaussian", opt = "rff", 
                      rho = rho, selection_method = "REML", verbose = TRUE)
lambda_rff_n5000d7 = FastKRR_rff$lambda


### n = 10000, d = 3
n = 10000
d = 3
X = matrix(runif(n*d,0, 1),n,d)
y = as.vector(sin(2*pi*rowMeans(X)^3) + rnorm(n, 0, 0.1))

data_n10000d3 = data.frame(X, y)

FastKRR_exact = fastkrr(data = data_n10000d3, response = "y", kernel = "gaussian", opt = "exact", 
                        rho = rho, selection_method = "REML", verbose = TRUE)
lambda_exact_n10000d3 = FastKRR_exact$lambda


FastKRR_nystrom = fastkrr(data = data_n10000d3, response = "y", kernel = "gaussian", opt = "nystrom", 
                          rho = rho, selection_method = "REML", verbose = TRUE)
lambda_nystrom_n10000d3 = FastKRR_nystrom$lambda


FastKRR_pivoted = fastkrr(data = data_n10000d3, response = "y", kernel = "gaussian", opt = "pivoted", 
                          rho = rho, selection_method = "REML", verbose = TRUE)
lambda_pivoted_n10000d3 = FastKRR_pivoted$lambda


FastKRR_rff = fastkrr(data = data_n10000d3, response = "y", kernel = "gaussian", opt = "rff", 
                      rho = rho, selection_method = "REML", verbose = TRUE)
lambda_rff_n10000d3 = FastKRR_rff$lambda


### n = 10000, d = 5
n = 10000
d = 5
X = matrix(runif(n*d,0, 1),n,d)
y = as.vector(sin(2*pi*rowMeans(X)^3) + rnorm(n, 0, 0.1))

data_n10000d5 = data.frame(X, y)

FastKRR_exact = fastkrr(data = data_n10000d5, response = "y", kernel = "gaussian", opt = "exact", 
                        rho = rho, selection_method = "REML", verbose = TRUE)
lambda_exact_n10000d5 = FastKRR_exact$lambda


FastKRR_nystrom = fastkrr(data = data_n10000d5, response = "y", kernel = "gaussian", opt = "nystrom", 
                          rho = rho, selection_method = "REML", verbose = TRUE)
lambda_nystrom_n10000d5 = FastKRR_nystrom$lambda


FastKRR_pivoted = fastkrr(data = data_n10000d5, response = "y", kernel = "gaussian", opt = "pivoted", 
                          rho = rho, selection_method = "REML", verbose = TRUE)
lambda_pivoted_n10000d5 = FastKRR_pivoted$lambda


FastKRR_rff = fastkrr(data = data_n10000d5, response = "y", kernel = "gaussian", opt = "rff", 
                      rho = rho, selection_method = "REML", verbose = TRUE)
lambda_rff_n10000d5 = FastKRR_rff$lambda


### n = 10000, d = 7
n = 10000
d = 7
X = matrix(runif(n*d,0, 1),n,d)
y = as.vector(sin(2*pi*rowMeans(X)^3) + rnorm(n, 0, 0.1))

data_n10000d7 = data.frame(X, y)

FastKRR_exact = fastkrr(data = data_n10000d7, response = "y", kernel = "gaussian", opt = "exact", 
                        rho = rho, selection_method = "REML", verbose = TRUE)
lambda_exact_n10000d7 = FastKRR_exact$lambda


FastKRR_nystrom = fastkrr(data = data_n10000d7, response = "y", kernel = "gaussian", opt = "nystrom", 
                          rho = rho, selection_method = "REML", verbose = TRUE)
lambda_nystrom_n10000d7 = FastKRR_nystrom$lambda


FastKRR_pivoted = fastkrr(data = data_n10000d7, response = "y", kernel = "gaussian", opt = "pivoted", 
                          rho = rho, selection_method = "REML", verbose = TRUE)
lambda_pivoted_n10000d7 = FastKRR_pivoted$lambda


FastKRR_rff = fastkrr(data = data_n10000d7, response = "y", kernel = "gaussian", opt = "rff", 
                      rho = rho, selection_method = "REML", verbose = TRUE)
lambda_rff_n10000d7 = FastKRR_rff$lambda



# time table
time_data = microbenchmark::microbenchmark(
  FastKRR_exact_n5000d3 = fastkrr(data = data_n5000d3, response = "y", kernel = "gaussian", opt = "exact", 
                                  n_threads = 4, rho = rho, lambda = lambda_exact_n5000d3, verbose = FALSE),
  
  FastKRR_nystrom_n5000d3 = fastkrr(data = data_n5000d3, response = "y", kernel = "gaussian", opt = "nystrom", 
                                    n_threads = 4, rho = rho, lambda = lambda_nystrom_n5000d3, verbose = FALSE),
  
  FastKRR_pivoted_n5000d3 = fastkrr(data = data_n5000d3, response = "y", kernel = "gaussian", opt = "pivoted",
                                    n_threads = 4, rho = rho, lambda = lambda_pivoted_n5000d3, verbose = FALSE),
  
  FastKRR_rff_n5000d3 = fastkrr(data = data_n5000d3, response = "y", kernel = "gaussian", opt = "rff", 
                                n_threads = 4, rho = rho, lambda = lambda_rff_n5000d3, verbose = FALSE),
  
  
  
  
  FastKRR_exact_n5000d5 = fastkrr(data = data_n5000d5, response = "y", kernel = "gaussian", opt = "exact", 
                                  n_threads = 4, rho = rho, lambda = lambda_exact_n5000d5, verbose = FALSE),
  
  FastKRR_nystrom_n5000d5 = fastkrr(data = data_n5000d5, response = "y", kernel = "gaussian", opt = "nystrom", 
                                    n_threads = 4, rho = rho, lambda = lambda_nystrom_n5000d5, verbose = FALSE),
  
  FastKRR_pivoted_n5000d5 = fastkrr(data = data_n5000d5, response = "y", kernel = "gaussian", opt = "pivoted",
                                    rho = rho, lambda = lambda_pivoted_n5000d5, verbose = FALSE),
  
  FastKRR_rff_n5000d5 = fastkrr(data = data_n5000d5, response = "y", kernel = "gaussian", opt = "rff", 
                                n_threads = 4, rho = rho, lambda = lambda_rff_n5000d5, verbose = FALSE),
  
  
  
  
  FastKRR_exact_n5000d7 = fastkrr(data = data_n5000d7, response = "y", kernel = "gaussian", opt = "exact", 
                                  n_threads = 4, rho = rho, lambda = lambda_exact_n5000d7, verbose = FALSE),
  
  FastKRR_nystrom_n5000d7 = fastkrr(data = data_n5000d7, response = "y", kernel = "gaussian", opt = "nystrom", 
                                    n_threads = 4, rho = rho, lambda = lambda_nystrom_n5000d7, verbose = FALSE),
  
  FastKRR_pivoted_n5000d7 = fastkrr(data = data_n5000d7, response = "y", kernel = "gaussian", opt = "pivoted",
                                    n_threads = 4, rho = rho, lambda = lambda_pivoted_n5000d7, verbose = FALSE),
  
  FastKRR_rff_n5000d7 = fastkrr(data = data_n5000d7, response = "y", kernel = "gaussian", opt = "rff", 
                                n_threads = 4, rho = rho, lambda = lambda_rff_n5000d7, verbose = FALSE),
  
  
  
  
  FastKRR_exact_n10000d3 = fastkrr(data = data_n10000d3, response = "y", kernel = "gaussian", opt = "exact", 
                                   n_threads = 4, rho = rho, lambda = lambda_exact_n10000d3, verbose = FALSE),
  
  FastKRR_nystrom_n10000d3 = fastkrr(data = data_n10000d3, response = "y", kernel = "gaussian", opt = "nystrom", 
                                     n_threads = 4, rho = rho, lambda = lambda_nystrom_n10000d3, verbose = FALSE),
  
  FastKRR_pivoted_n10000d3 = fastkrr(data = data_n10000d3, response = "y", kernel = "gaussian", opt = "pivoted",
                                     n_threads = 4, rho = rho, lambda = lambda_pivoted_n10000d3, verbose = FALSE),
  
  FastKRR_rff_n10000d3 = fastkrr(data = data_n10000d3, response = "y", kernel = "gaussian", opt = "rff", 
                                 n_threads = 4, rho = rho, lambda = lambda_rff_n10000d3, verbose = FALSE),
  
  
  
  FastKRR_exact_n10000d5 = fastkrr(data = data_n10000d5, response = "y", kernel = "gaussian", opt = "exact", 
                                   n_threads = 4, rho = rho, lambda = lambda_exact_n10000d5, verbose = FALSE),
  
  FastKRR_nystrom_n10000d5 = fastkrr(data = data_n10000d5, response = "y", kernel = "gaussian", opt = "nystrom", 
                                     n_threads = 4, rho = rho, lambda = lambda_nystrom_n10000d5, verbose = FALSE),
  
  FastKRR_pivoted_n10000d5 = fastkrr(data = data_n10000d5, response = "y", kernel = "gaussian", opt = "pivoted",
                                     n_threads = 4, rho = rho, lambda = lambda_pivoted_n10000d5, verbose = FALSE),
  
  FastKRR_rff_n10000d5 = fastkrr(data = data_n10000d5, response = "y", kernel = "gaussian", opt = "rff", 
                                 n_threads = 4, rho = rho, lambda = lambda_rff_n10000d5, verbose = FALSE),
  
  
  
  
  FastKRR_exact_n10000d7 = fastkrr(data = data_n10000d7, response = "y", kernel = "gaussian", opt = "exact", 
                                   n_threads = 4, rho = rho, lambda = lambda_exact_n10000d7, verbose = FALSE),
  
  FastKRR_nystrom_n10000d7 = fastkrr(data = data_n10000d7, response = "y", kernel = "gaussian", opt = "nystrom", 
                                     n_threads = 4, rho = rho, lambda = lambda_nystrom_n10000d7, verbose = FALSE),
  
  FastKRR_pivoted_n10000d7 = fastkrr(data = data_n10000d7, response = "y", kernel = "gaussian", opt = "pivoted",
                                     n_threads = 4, rho = rho, lambda = lambda_pivoted_n10000d7, verbose = FALSE),
  
  FastKRR_rff_n10000d7 = fastkrr(data = data_n10000d7, response = "y", kernel = "gaussian", opt = "rff", 
                                 n_threads = 4, rho = rho, lambda = lambda_rff_n10000d7, verbose = FALSE),
  times = 100,
  unit = "ms"
)

save(time_data, file = "output/table3_compare_time_opt(all).Rdata")

time_tbl = time_data %>%
  as.data.frame() %>%
  mutate(
    expr = as.character(expr),
    Method = sub("_n[0-9]+d[0-9]+$", "", expr),
    n = as.integer(sub(".*_n([0-9]+)d[0-9]+$", "\\1", expr)),
    d = as.integer(sub(".*d([0-9]+)$", "\\1", expr)),
    time_sec = time / 1e9
  ) %>%
  group_by(Method, n, d) %>%
  summarise(
    Mean = mean(time_sec, na.rm = TRUE),
    SE = sd(time_sec, na.rm = TRUE) /
      sqrt(sum(!is.na(time_sec))),
    Median = median(time_sec, na.rm = TRUE),
    .groups = "drop"
  )

# prediction
## n 5000, d 3
compare_packages_exact_n5000d3 = list()

for(rep in 1:100){
  n = 5000
  d = 3
  X = matrix(runif(n*d,0, 1),n,d)
  y = as.vector(sin(2*pi*rowMeans(X)^3) + rnorm(n, 0, 0.1))
  
  data = data.frame(X, y)
  
  n_z = ceiling(1.2 * n)
  
  z = matrix(seq(0, 1, len = n_z * d),n_z, d)
  y_z = as.vector(sin(2*pi*rowMeans(z)^3))
  
  FastKRR_exact = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "exact", 
                          rho = rho, lambda = lambda_exact_n5000d3, verbose = FALSE)
  
  FastKRR_nystrom = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "nystrom", 
                            rho = rho, lambda = lambda_nystrom_n5000d3, verbose = FALSE)
  
  FastKRR_pivoted = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "pivoted",
                            rho = rho, lambda = lambda_pivoted_n5000d3, verbose = FALSE)
  
  FastKRR_rff = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "rff", 
                        rho = rho, lambda = lambda_rff_n5000d3, verbose = FALSE)
  
  
  compare_packages_exact_n5000d3$exact[[rep]] =  as.vector(y_z - predict(object = FastKRR_exact, newdata = z))
  compare_packages_exact_n5000d3$nystrom[[rep]] =  as.vector(y_z - predict(object = FastKRR_nystrom, newdata = z))
  compare_packages_exact_n5000d3$pivoted[[rep]] =  as.vector(y_z - predict(object = FastKRR_pivoted, newdata = z))
  compare_packages_exact_n5000d3$rff[[rep]] =  as.vector(y_z - predict(object = FastKRR_rff, newdata = z))
}

save(compare_packages_exact_n5000d3, file = "output/table3_compare_prediction (n5000d3).Rdata")


## n5000, d5
compare_packages_exact_n5000d5 = list()

for(rep in 1:100){
  n = 5000
  d = 5
  X = matrix(runif(n*d,0, 1),n,d)
  y = as.vector(sin(2*pi*rowMeans(X)^3) + rnorm(n, 0, 0.1))
  
  data = data.frame(X, y)
  
  n_z = ceiling(1.2 * n)
  
  z = matrix(seq(0, 1, len = n_z * d),n_z, d)
  y_z = as.vector(sin(2*pi*rowMeans(z)^3))
  
  FastKRR_exact = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "exact", 
                          rho = rho, lambda = lambda_exact_n5000d5, verbose = FALSE)
  
  FastKRR_nystrom = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "nystrom", 
                            rho = rho, lambda = lambda_nystrom_n5000d5, verbose = FALSE)
  
  FastKRR_pivoted = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "pivoted",
                            rho = rho, lambda = lambda_pivoted_n5000d5, verbose = FALSE)
  
  FastKRR_rff = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "rff", 
                        rho = rho, lambda = lambda_rff_n5000d5, verbose = FALSE)
  
  
  
  
  compare_packages_exact_n5000d5$exact[[rep]] =  as.vector(y_z - predict(object = FastKRR_exact, newdata = z))
  compare_packages_exact_n5000d5$nystrom[[rep]] =  as.vector(y_z - predict(object = FastKRR_nystrom, newdata = z))
  compare_packages_exact_n5000d5$pivoted[[rep]] =  as.vector(y_z - predict(object = FastKRR_pivoted, newdata = z))
  compare_packages_exact_n5000d5$rff[[rep]] =  as.vector(y_z - predict(object = FastKRR_rff, newdata = z))
}

save(compare_packages_exact_n5000d5, file = "output/table3_compare_prediction (n5000d5).Rdata")


## n5000, d7
compare_packages_exact_n5000d7 = list()

for(rep in 1:100){
  n = 5000
  d = 7
  X = matrix(runif(n*d,0, 1),n,d)
  y = as.vector(sin(2*pi*rowMeans(X)^3) + rnorm(n, 0, 0.1))
  
  data = data.frame(X, y)
  
  n_z = ceiling(1.2 * n)
  
  z = matrix(seq(0, 1, len = n_z * d),n_z, d)
  y_z = as.vector(sin(2*pi*rowMeans(z)^3))
  
  FastKRR_exact = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "exact", 
                          rho = rho, lambda = lambda_exact_n5000d7, verbose = FALSE)
  
  FastKRR_nystrom = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "nystrom", 
                            rho = rho, lambda = lambda_nystrom_n5000d7, verbose = FALSE)
  
  FastKRR_pivoted = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "pivoted",
                            rho = rho, lambda = lambda_pivoted_n5000d7, verbose = FALSE)
  
  FastKRR_rff = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "rff", 
                        rho = rho, lambda = lambda_rff_n5000d7, verbose = FALSE)
  
  
  
  
  compare_packages_exact_n5000d7$exact[[rep]] =  as.vector(y_z - predict(object = FastKRR_exact, newdata = z))
  compare_packages_exact_n5000d7$nystrom[[rep]] =  as.vector(y_z - predict(object = FastKRR_nystrom, newdata = z))
  compare_packages_exact_n5000d7$pivoted[[rep]] =  as.vector(y_z - predict(object = FastKRR_pivoted, newdata = z))
  compare_packages_exact_n5000d7$rff[[rep]] =  as.vector(y_z - predict(object = FastKRR_rff, newdata = z))
}

save(compare_packages_exact_n5000d7, file = "output/table3_compare_prediction (n5000d7).Rdata")


## n10000, d3
compare_packages_exact_n10000d3 = list()

for(rep in 1:100){
  n = 10000
  d = 3
  X = matrix(runif(n*d,0, 1),n,d)
  y = as.vector(sin(2*pi*rowMeans(X)^3) + rnorm(n, 0, 0.1))
  
  data = data.frame(X, y)
  
  n_z = ceiling(1.2 * n)
  
  z = matrix(seq(0, 1, len = n_z * d),n_z, d)
  y_z = as.vector(sin(2*pi*rowMeans(z)^3))
  
  FastKRR_exact = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "exact", 
                          rho = rho, lambda = lambda_exact_n10000d3, verbose = FALSE)
  
  FastKRR_nystrom = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "nystrom", 
                            rho = rho, lambda = lambda_nystrom_n10000d3, verbose = FALSE)
  
  FastKRR_pivoted = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "pivoted",
                            rho = rho, lambda = lambda_pivoted_n10000d3, verbose = FALSE)
  
  FastKRR_rff = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "rff", 
                        rho = rho, lambda = lambda_rff_n10000d3, verbose = FALSE)
  
  
  compare_packages_exact_n10000d3$exact[[rep]] =  as.vector(y_z - predict(object = FastKRR_exact, newdata = z))
  compare_packages_exact_n10000d3$nystrom[[rep]] =  as.vector(y_z - predict(object = FastKRR_nystrom, newdata = z))
  compare_packages_exact_n10000d3$pivoted[[rep]] =  as.vector(y_z - predict(object = FastKRR_pivoted, newdata = z))
  compare_packages_exact_n10000d3$rff[[rep]] =  as.vector(y_z - predict(object = FastKRR_rff, newdata = z))
}

save(compare_packages_exact_n10000d3, file = "output/table3_compare_prediction (n10000d3).Rdata")


## n10000, d5
compare_packages_exact_n10000d5 = list()

for(rep in 1:100){
  n = 10000
  d = 5
  X = matrix(runif(n*d,0, 1),n,d)
  y = as.vector(sin(2*pi*rowMeans(X)^3) + rnorm(n, 0, 0.1))
  
  data = data.frame(X, y)
  
  n_z = ceiling(1.2 * n)
  
  z = matrix(seq(0, 1, len = n_z * d),n_z, d)
  y_z = as.vector(sin(2*pi*rowMeans(z)^3))
  
  FastKRR_exact = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "exact", 
                          rho = rho, lambda = lambda_exact_n10000d5, verbose = FALSE)
  
  FastKRR_nystrom = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "nystrom", 
                            rho = rho, lambda = lambda_nystrom_n10000d5, verbose = FALSE)
  
  FastKRR_pivoted = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "pivoted",
                            rho = rho, lambda = lambda_pivoted_n10000d5, verbose = FALSE)
  
  FastKRR_rff = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "rff", 
                        rho = rho, lambda = lambda_rff_n10000d5, verbose = FALSE)
  
  
  
  
  compare_packages_exact_n10000d5$exact[[rep]] =  as.vector(y_z - predict(object = FastKRR_exact, newdata = z))
  compare_packages_exact_n10000d5$nystrom[[rep]] =  as.vector(y_z - predict(object = FastKRR_nystrom, newdata = z))
  compare_packages_exact_n10000d5$pivoted[[rep]] =  as.vector(y_z - predict(object = FastKRR_pivoted, newdata = z))
  compare_packages_exact_n10000d5$rff[[rep]] =  as.vector(y_z - predict(object = FastKRR_rff, newdata = z))
}

save(compare_packages_exact_n10000d5, file = "output/table3_compare_prediction (n10000d5).Rdata")


## n10000, d7
compare_packages_exact_n10000d7 = list()

for(rep in 1:100){
  n = 10000
  d = 7
  X = matrix(runif(n*d,0, 1),n,d)
  y = as.vector(sin(2*pi*rowMeans(X)^3) + rnorm(n, 0, 0.1))
  
  data = data.frame(X, y)
  
  n_z = ceiling(1.2 * n)
  
  z = matrix(seq(0, 1, len = n_z * d),n_z, d)
  y_z = as.vector(sin(2*pi*rowMeans(z)^3))
  
  FastKRR_exact = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "exact", 
                          rho = rho, lambda = lambda_exact_n10000d7, verbose = FALSE)
  
  FastKRR_nystrom = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "nystrom", 
                            rho = rho, lambda = lambda_nystrom_n10000d7, verbose = FALSE)
  
  FastKRR_pivoted = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "pivoted",
                            rho = rho, lambda = lambda_pivoted_n10000d7, verbose = FALSE)
  
  FastKRR_rff = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "rff", 
                        rho = rho, lambda = lambda_rff_n10000d7, verbose = FALSE)
  
  
  
  
  compare_packages_exact_n10000d7$exact[[rep]] =  as.vector(y_z - predict(object = FastKRR_exact, newdata = z))
  compare_packages_exact_n10000d7$nystrom[[rep]] =  as.vector(y_z - predict(object = FastKRR_nystrom, newdata = z))
  compare_packages_exact_n10000d7$pivoted[[rep]] =  as.vector(y_z - predict(object = FastKRR_pivoted, newdata = z))
  compare_packages_exact_n10000d7$rff[[rep]] =  as.vector(y_z - predict(object = FastKRR_rff, newdata = z))
}

save(compare_packages_exact_n10000d7, file = "output/table3_compare_prediction (n10000d7).Rdata")

summary_metric_n5000d3 = do.call(rbind, lapply(compare_packages_exact_n5000d3, calc_metric))

summary_metric_n5000d5 = do.call(rbind, lapply(compare_packages_exact_n5000d5, calc_metric))

summary_metric_n5000d7 = do.call(rbind, lapply(compare_packages_exact_n5000d7, calc_metric))

summary_metric_n10000d3 = do.call(rbind, lapply(compare_packages_exact_n10000d3, calc_metric))

summary_metric_n10000d5 = do.call(rbind, lapply(compare_packages_exact_n10000d5, calc_metric))

summary_metric_n10000d7 = do.call(rbind, lapply(compare_packages_exact_n10000d7, calc_metric))

metric_tbl = bind_rows(
  summary_metric_n5000d3 %>%
    rownames_to_column("method_short") %>%
    mutate(n = 5000, d = 3),
  
  summary_metric_n5000d5 %>%
    rownames_to_column("method_short") %>%
    mutate(n = 5000, d = 5),
  
  summary_metric_n5000d7 %>%
    rownames_to_column("method_short") %>%
    mutate(n = 5000, d = 7),
  
  summary_metric_n10000d3 %>%
    rownames_to_column("method_short") %>%
    mutate(n = 10000, d = 3),
  
  summary_metric_n10000d5 %>%
    rownames_to_column("method_short") %>%
    mutate(n = 10000, d = 5),
  
  summary_metric_n10000d7 %>%
    rownames_to_column("method_short") %>%
    mutate(n = 10000, d = 7)
) %>%
  mutate(
    Method = recode(
      method_short,
      exact   = "FastKRR_exact",
      nystrom = "FastKRR_nystrom",
      pivoted = "FastKRR_pivoted",
      rff     = "FastKRR_rff"
    )
  ) %>%
  select(
    Method,
    n,
    d,
    MSE_mean,
    MSE_se,
    MAE_mean,
    MAE_se
  )

##total table 
final_tbl = time_tbl %>%
  left_join(metric_tbl, by = c("Method", "n", "d")) %>%
  mutate(
    Method = factor(
      Method,
      levels = c(
        "FastKRR_exact",
        "FastKRR_nystrom",
        "FastKRR_pivoted",
        "FastKRR_rff"
      )
    )
  ) %>%
  arrange(Method, n, d)


final_tbl = time_tbl %>%
  left_join(metric_tbl, by = c("Method", "n", "d")) %>%
  mutate(Method = factor(Method, levels = c(
    "FastKRR_exact",
    "FastKRR_nystrom",
    "FastKRR_pivoted",
    "FastKRR_rff"))) %>%
  arrange(Method, n, d)

final_tbl_print = final_tbl %>%
  arrange(Method, n, d) %>%
  mutate(
    Method = recode(
      as.character(Method),
      "FastKRR_exact"    = "exact",
      "FastKRR_nystrom" = "nystrom",
      "FastKRR_pivoted" = "pivoted",
      "FastKRR_rff"     = "rff"
    )) %>%
  group_by(Method) %>%
  mutate(
    Method_print = ifelse(row_number() == 1, paste0("\\textbf{", Method, "}"), "~"),
    
    n_print = ifelse(!duplicated(n), as.character(n), "~"),
    
    `Mean (SE)` = sprintf("%.4f (%.4f)", Mean, SE),
    
    Median = sprintf("%.4f", Median),
    
    `MSE (SE)` = sprintf("%.4f (%.4f)", MSE_mean, MSE_se),
    
    `MAE (SE)` = sprintf("%.4f (%.4f)", MAE_mean, MAE_se)) %>%
  ungroup() %>%
  select(
    Method = Method_print,
    n = n_print,
    d,
    `Mean (SE)`,
    Median,
    `MSE (SE)`,
    `MAE (SE)`
  )


## 6.3 Efficiency of the approximation methods in FastKR
#m = 111
#data generate
n = 10000
d = 3
X = matrix(runif(n*d, 0, 1), n, d)
colnames(X) = c("X1", "X2", "X3")
X1 = X[, 1]; X2 = X[, 2]; X3 = X[, 3]
y = as.vector(sin(2*pi*rowMeans(X)^3) + rnorm(n, 0, 0.1))
data = data.frame(X, y = y)

m = 111
# ---- KRLS: LOOCV로 튜닝 ----
model_krls_tune = krls(X = X, y = y, whichkernel = "gaussian",
                       sigma = rho, derivative = FALSE,
                       binary = FALSE, vcov = FALSE,
                       approx = "nystrom", nystrom_m = m)

lambda_krls_n10000d3 = model_krls_tune$lambda

# ---- gKRLS: REML로 튜닝 ----
model_gkrls_tune = gam(y ~ s(X1, X2, X3, bs = "gKRLS",
                             xt = gKRLS(
                               sketch_method = "subsampling",
                               sketch_multiplier = NULL,
                               sketch_size_raw = m,
                               remove_instability = TRUE
                             )),
                       data = data,
                       method = "REML")

sp_gkrls_n10000d3 = model_gkrls_tune$sp


# ---- FastKRR: CVST(REML)로 튜닝 ----
FastKRR_nystrom = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "nystrom",
                          rho = rho, selection_method = "REML", verbose = TRUE)

lambda_nystrom_n10000d3 = FastKRR_nystrom$lambda


FastKRR_pivoted = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "pivoted",
                          rho = rho, selection_method = "REML", verbose = TRUE)

lambda_pivoted_n10000d3 = FastKRR_pivoted$lambda


FastKRR_rff = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "rff",
                      rho = rho, selection_method = "REML", verbose = TRUE)

lambda_rff_n10000d3 = FastKRR_rff$lambda

##time
## 시간
krr_pkg_compare_n10000d3 = microbenchmark::microbenchmark(
  KRLS = krls(X = X, y = y, whichkernel = "gaussian",
              sigma = rho, derivative = FALSE, lambda = lambda_krls_n10000d3,
              binary = FALSE, vcov = FALSE, approx = "nystrom", nystrom_m = m),
  
  gKRLS = gam(y ~ s(X1, X2, X3, bs = "gKRLS",
                    xt = gKRLS(
                      sketch_method = "subsampling",
                      sketch_size_raw = 111,
                      sketch_multiplier = NULL,
                      remove_instability = TRUE
                    )),
              data = data,
              sp = sp_gkrls_n10000d3),
  
  FastKRR_nystrom = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "nystrom",
                            rho = rho, lambda = lambda_nystrom_n10000d3, verbose = FALSE, m = 111),
  
  FastKRR_pivoted = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "pivoted",
                            rho = rho, lambda = lambda_pivoted_n10000d3, verbose = FALSE, m = 111),
  
  FastKRR_rff = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "rff",
                        rho = rho, lambda = lambda_rff_n10000d3, verbose = FALSE, m = 111),
  times = 100,
  unit = "ms"
)

save(krr_pkg_compare_n10000d3, file = "output/table4_compare_time_pkg(n10000d3_m111).Rdata")

krr_pkg_compare_n10000d3m111 = krr_pkg_compare_n10000d3

## pred

compare_packages_approx_n10000d3 = list()

for (rep in 1:100) {
  n = 10000
  d = 3
  X = matrix(runif(n*d, 0, 1), n, d)
  colnames(X) = c("X1", "X2", "X3")
  X1 = X[, 1]; X2 = X[, 2]; X3 = X[, 3]
  y = as.vector(sin(2*pi*rowMeans(X)^3) + rnorm(n, 0, 0.1))
  data = data.frame(X, y = y)
  
  n_z = ceiling(1.2 * n)
  z = matrix(runif(n_z*d, 0, 1), n_z, d)
  colnames(z) = c("X1", "X2", "X3")
  y_z = as.vector(sin(2*pi*rowMeans(z)^3))
  z_frame = as.data.frame(z)
  
  # ---- KRLS (lambda 고정) ----
  model_krls = krls(X = X, y = y, whichkernel = "gaussian",
                    sigma = rho, derivative = FALSE, lambda = lambda_krls_n10000d3,
                    binary = FALSE, vcov = FALSE, approx = "nystrom", nystrom_m = 111)
  pred_krls = predict(object = model_krls, newdata = as.matrix(z))$fit
  
  # ---- gKRLS (sp 고정) ----
  model_gkrls = gam(y ~ s(X1, X2, X3, bs = "gKRLS",
                          xt = gKRLS(
                            sketch_method = "subsampling",
                            sketch_size_raw = 111,
                            sketch_multiplier = NULL,
                            remove_instability = TRUE
                          )),
                    data = data,
                    sp = sp_gkrls_n10000d3)
  pred_gkrls = predict(object = model_gkrls, newdata = z_frame)
  
  # ---- FastKRR (lambda 고정, fastcv 없이) ----
  FastKRR_nystrom = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "nystrom",
                            rho = rho, lambda = lambda_nystrom_n10000d3, verbose = FALSE, m = 111)
  
  FastKRR_pivoted = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "pivoted",
                            rho = rho, lambda = lambda_pivoted_n10000d3, verbose = FALSE, m = 111)
  
  FastKRR_rff = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "rff",
                        rho = rho, lambda = lambda_rff_n10000d3, verbose = FALSE, m = 111)
  
  # ---- 저장 ----
  compare_packages_approx_n10000d3$KRLS[[rep]]           = as.vector(y_z - pred_krls)
  compare_packages_approx_n10000d3$gKRLS[[rep]]           = as.vector(y_z - pred_gkrls)
  compare_packages_approx_n10000d3$FastKRR_nystrom[[rep]] = as.vector(y_z - predict(object = FastKRR_nystrom, newdata = as.matrix(z_frame)))
  compare_packages_approx_n10000d3$FastKRR_pivoted[[rep]] = as.vector(y_z - predict(object = FastKRR_pivoted, newdata = as.matrix(z_frame)))
  compare_packages_approx_n10000d3$FastKRR_rff[[rep]]     = as.vector(y_z - predict(object = FastKRR_rff, newdata = as.matrix(z_frame)))
}

save(compare_packages_approx_n10000d3, file = "output/table4_compare_prediction(n10000d3_m111).Rdata")

compare_packages_approx_n10000d3m111 = compare_packages_approx_n10000d3

# m = 208
n = 10000
d = 3
X = matrix(runif(n*d, 0, 1), n, d)
colnames(X) = c("X1", "X2", "X3")
X1 = X[, 1]; X2 = X[, 2]; X3 = X[, 3]
y = as.vector(sin(2*pi*rowMeans(X)^3) + rnorm(n, 0, 0.1))
data = data.frame(X, y = y)

m = ceiling(n^(1/2) * log(d + 5))
m

# ---- KRLS: LOOCV로 튜닝 ----
model_krls_tune = krls(X = X, y = y, whichkernel = "gaussian",
                       sigma = rho, derivative = FALSE,
                       binary = FALSE, vcov = FALSE,
                       approx = "nystrom", nystrom_m = m)

lambda_krls_n10000d3 = model_krls_tune$lambda


# ---- gKRLS: REML로 튜닝 ----
model_gkrls_tune = gam(y ~ s(X1, X2, X3, bs = "gKRLS",
                             xt = gKRLS(
                               sketch_method = "subsampling",
                               sketch_multiplier = NULL,
                               sketch_size_raw = m,
                               remove_instability = TRUE
                             )),
                       data = data,
                       method = "REML")

sp_gkrls_n10000d3 = model_gkrls_tune$sp


# ---- FastKRR: CVST(REML)로 튜닝 ----
FastKRR_nystrom = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "nystrom",
                          rho = rho, selection_method = "REML", verbose = TRUE)

lambda_nystrom_n10000d3 = FastKRR_nystrom$lambda


FastKRR_pivoted = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "pivoted",
                          rho = rho, selection_method = "REML", verbose = TRUE)

lambda_pivoted_n10000d3 = FastKRR_pivoted$lambda


FastKRR_rff = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "rff",
                      rho = rho, selection_method = "REML", verbose = TRUE)

lambda_rff_n10000d3 = FastKRR_rff$lambda

## time
krr_pkg_compare_n10000d3 = microbenchmark::microbenchmark(
  KRLS = krls(X = X, y = y, whichkernel = "gaussian",
              sigma = rho, derivative = FALSE, lambda = lambda_krls_n10000d3,
              binary = FALSE, vcov = FALSE, approx = "nystrom", nystrom_m = 208),
  
  gKRLS = gam(y ~ s(X1, X2, X3, bs = "gKRLS",
                    xt = gKRLS(
                      sketch_method = "subsampling",
                      sketch_size_raw = 208,
                      sketch_multiplier = NULL,
                      remove_instability = TRUE
                    )),
              data = data,
              sp = sp_gkrls_n10000d3),
  
  FastKRR_nystrom = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "nystrom",
                            rho = rho, lambda = lambda_nystrom_n10000d3, verbose = FALSE, m = 208),
  
  FastKRR_pivoted = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "pivoted",
                            rho = rho, lambda = lambda_pivoted_n10000d3, verbose = FALSE, m = 208),
  
  FastKRR_rff = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "rff",
                        rho = rho, lambda = lambda_rff_n10000d3, verbose = FALSE, m = 208),
  times = 100,
  unit = "ms"
)
save(krr_pkg_compare_n10000d3, file = "output/table4_compare_time_pkg(n10000d3_m208).Rdata")

krr_pkg_compare_n10000d3m208 = krr_pkg_compare_n10000d3

## pred

compare_packages_approx_n10000d3 = list()

for (rep in 1:100) {
  n = 10000
  d = 3
  X = matrix(runif(n*d, 0, 1), n, d)
  colnames(X) = c("X1", "X2", "X3")
  X1 = X[, 1]; X2 = X[, 2]; X3 = X[, 3]
  y = as.vector(sin(2*pi*rowMeans(X)^3) + rnorm(n, 0, 0.1))
  data = data.frame(X, y = y)
  
  n_z = ceiling(1.2 * n)
  z = matrix(runif(n_z*d, 0, 1), n_z, d)
  colnames(z) = c("X1", "X2", "X3")
  y_z = as.vector(sin(2*pi*rowMeans(z)^3))
  z_frame = as.data.frame(z)
  
  # ---- KRLS (lambda 고정) ----
  model_krls = krls(X = X, y = y, whichkernel = "gaussian",
                    sigma = rho, derivative = FALSE, lambda = lambda_krls_n10000d3,
                    binary = FALSE, vcov = FALSE, approx = "nystrom", nystrom_m = 208)
  pred_krls = predict(object = model_krls, newdata = as.matrix(z))$fit
  
  # ---- gKRLS (sp 고정) ----
  model_gkrls = gam(y ~ s(X1, X2, X3, bs = "gKRLS",
                          xt = gKRLS(
                            sketch_method = "subsampling",
                            sketch_size_raw = 208,
                            sketch_multiplier = NULL,
                            remove_instability = TRUE
                          )),
                    data = data,
                    sp = sp_gkrls_n10000d3)
  pred_gkrls = predict(object = model_gkrls, newdata = z_frame)
  
  # ---- FastKRR (lambda 고정, fastcv 없이) ----
  FastKRR_nystrom = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "nystrom",
                            rho = rho, lambda = lambda_nystrom_n10000d3, verbose = FALSE, m = 208)
  
  FastKRR_pivoted = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "pivoted",
                            rho = rho, lambda = lambda_pivoted_n10000d3, verbose = FALSE, m = 208)
  
  FastKRR_rff = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "rff",
                        rho = rho, lambda = lambda_rff_n10000d3, verbose = FALSE, m = 208)
  
  # ---- 저장 ----
  compare_packages_approx_n10000d3$KRLS[[rep]]           = as.vector(y_z - pred_krls)
  compare_packages_approx_n10000d3$gKRLS[[rep]]           = as.vector(y_z - pred_gkrls)
  compare_packages_approx_n10000d3$FastKRR_nystrom[[rep]] = as.vector(y_z - predict(object = FastKRR_nystrom, newdata = as.matrix(z_frame)))
  compare_packages_approx_n10000d3$FastKRR_pivoted[[rep]] = as.vector(y_z - predict(object = FastKRR_pivoted, newdata = as.matrix(z_frame)))
  compare_packages_approx_n10000d3$FastKRR_rff[[rep]]     = as.vector(y_z - predict(object = FastKRR_rff, newdata = as.matrix(z_frame)))
}

save(compare_packages_approx_n10000d3, file = "output/table4_compare_prediction(n10000d3_m208).Rdata")

compare_packages_approx_n10000d3m208 = compare_packages_approx_n10000d3


names(compare_packages_approx_n10000d3m111) = 
  c("KRLS-nystrom", "gKRLS-sketched", "FastKRR-nystrom", "FastKRR-pivoted", "FastKRR-rff")
krr_pkg_compare_n10000d3m111$expr = factor(
  krr_pkg_compare_n10000d3m111$expr,
  levels = c("gKRLS", "FastKRR_nystrom", "FastKRR_rff", "FastKRR_pivoted", "KRLS"),
  labels = c("gKRLS-sketched", "FastKRR-nystrom", "FastKRR-rff", "FastKRR-pivoted", "KRLS-nystrom")
)

names(compare_packages_approx_n10000d3m208) = 
  c("KRLS-nystrom", "gKRLS-sketched", "FastKRR-nystrom", "FastKRR-pivoted", "FastKRR-rff")

krr_pkg_compare_n10000d3m208$expr = factor(
  krr_pkg_compare_n10000d3m208$expr,
  levels = c("gKRLS", "FastKRR_nystrom", "FastKRR_rff", "FastKRR_pivoted", "KRLS"),
  labels = c("gKRLS-sketched", "FastKRR-nystrom", "FastKRR-rff", "FastKRR-pivoted", "KRLS-nystrom")
)


table_m111 = summarize_package_comparison(krr_pkg_compare_n10000d3m111, 
                                          compare_packages_approx_n10000d3m111, 111)

table_m208 = summarize_package_comparison(krr_pkg_compare_n10000d3m208, 
                                          compare_packages_approx_n10000d3m208, 208)


method_order = c(
  "KRLS-nystrom",
  "gKRLS-sketched",
  "FastKRR-nystrom",
  "FastKRR-pivoted",
  "FastKRR-rff"
)


table_final = bind_rows(
  table_m111,
  table_m208
) %>%
  mutate(
    Method = factor(
      Method,
      levels = method_order
    )
  ) %>%
  arrange(m, Method) %>%
  mutate(
    Method = as.character(Method),
    m_display = ifelse(
      duplicated(m),
      "",
      as.character(m)
    )
  ) %>%
  select(
    `m` = m_display,
    Method,
    `Mean (SE)`,
    Median,
    `MSE (SE)`,
    `MAE (SE)`
  )


## 6.4 Model selection efficiency
#simualtion---------------------------------------------------------------------
set.seed(1)
n = 2000 
d = 3
rho = 1
lambda = 10^seq(-11, -7, length.out = 100)
lambda_reml = c(10^(-11), 10^(-7))

X = matrix(runif(n*d, 0, 1), nrow = n, ncol = d)
y = as.vector(sin(2*pi*rowMeans(X)^3) + rnorm(n, 0, 0.1))
data_n2000d3 = data.frame(X, y)

# time
time_data_cv_n2000d3 = microbenchmark::microbenchmark(
  FastKRR_exactcv_n2000d3 = fastkrr(data = data_n2000d3, response = "y", kernel = "gaussian", opt = "exact", 
                                    n_threads = 4, rho = rho, lambda = lambda, selection_method = "exactCV", verbose = FALSE),
  
  FastKRR_fastcv_n2000d3 = fastkrr(data = data_n2000d3, response = "y", kernel = "gaussian", opt = "exact", 
                                   n_threads = 4, rho = rho, lambda = lambda, selection_method = "fastCV", verbose = FALSE),
  
  FastKRR_reml_n2000d3 = fastkrr(data = data_n2000d3, response = "y", kernel = "gaussian", opt = "exact", 
                                 n_threads = 4, rho = rho, lambda = lambda_reml, selection_method = "REML", verbose = FALSE),
  
  times = 100,
  unit = "ms"
)

save(time_data_cv_n2000d3, file = "output/table5_compare_time_cv.Rdata")

# pred
compare_cv_pred_n2000d3 = list()

for(rep in 1:100){
  n = 2000
  d = 3
  X = matrix(runif(n*d,0, 1),n,d)
  y = as.vector(sin(2*pi*rowMeans(X)^3) + rnorm(n, 0, 0.1))
  
  data = data.frame(X, y)
  
  n_z = ceiling(1.2 * n)
  z = matrix(seq(0, 1, len = n_z * d),n_z, d)
  y_z = as.vector(sin(2*pi*rowMeans(z)^3))
  
  FastKRR_exactCV = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "exact",
                            n_threads = 4, rho = rho, lambda = lambda, selection_method = "exactCV", verbose = TRUE)
  
  FastKRR_fastCV = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "exact", 
                           n_threads = 4, rho = rho, lambda = lambda, selection_method = "fastCV", verbose = TRUE)
  
  FastKRR_reml = fastkrr(data = data, response = "y", kernel = "gaussian", opt = "exact",
                         n_threads = 4, rho = rho, lambda = lambda_reml, selection_method = "REML", verbose = TRUE)
  
  compare_cv_pred_n2000d3$exactCV[[rep]] =  as.vector(y_z - predict(object = FastKRR_exactCV, newdata = z))
  compare_cv_pred_n2000d3$fastCV[[rep]] =  as.vector(y_z - predict(object = FastKRR_fastCV, newdata = z))
  compare_cv_pred_n2000d3$reml[[rep]] =  as.vector(y_z - predict(object = FastKRR_reml, newdata = z))
  
}

save(compare_cv_pred_n2000d3, file = "output/table5_compare_cv_prediction (n2000d3).Rdata")


df = as.data.frame(time_data_cv_n2000d3)


df$expr = factor(
  df$expr,
  levels = c("FastKRR_exactcv_n2000d3", "FastKRR_fastcv_n2000d3", "FastKRR_reml_n2000d3"),
  labels = c("ExactCV", "FastCV", "REML")
)

summary_cv = df %>%
  group_by(expr) %>%
  summarize(
    Mean     = mean(time) / 1e9,
    SE       = sd(time) / sqrt(100) / 1e9,
    Median   = median(time) / 1e9,
    .groups  = "drop"
  )

colnames(summary_cv)[colnames(summary_cv) == "expr"] = "Option"

summary_cv_pred_n2000d3 = do.call(rbind, lapply(compare_cv_pred_n2000d3, calc_metric))
summary_cv_table5 = cbind(summary_cv, summary_cv_pred_n2000d3)

summary_cv_table5 = summary_cv_table5 %>%mutate(
  `Mean (SE)` = sprintf("%.4f (%.4f)", Mean, SE),
  
  Median = sprintf(
    "%.4f",
    Median
  ),
  
  `MSE (SE)` = sprintf(
    "%.4f (%.4f)",
    MSE_mean,
    MSE_se
  ),
  
  `MAE (SE)` = sprintf(
    "%.4f (%.4f)",
    MAE_mean,
    MAE_se
  )
) %>%
  select(
    Option,
    `Mean (SE)`,
    Median,
    `MSE (SE)`,
    `MAE (SE)`
  )

rownames(summary_cv_table5) = NULL


## 7.Integration with the tidymodels ecosystem
## 7.4 Data preparation
data(ames)
ames = ames %>% mutate(Sale_Price = log10(Sale_Price))

set.seed(1)
ames_split = initial_split(ames, prop = 0.80, strata = Sale_Price)
ames_train = training(ames_split) # dim (2342, 74)
ames_test  = testing(ames_split) # dim (588, 74)

## 7.5 Model specification and workflow
# Model specification: exact KRR with Gaussian kernel
krr_spec = krr_reg(kernel = "gaussian", opt = "exact",
                   m = 50, eps = 1e-6, n_threads = 4,
                   rho = 1, penalty = tune()) %>%
  set_engine("fastkrr") %>%
  set_mode("regression")

# Recipe and workflow
rec = recipe(Sale_Price ~ Longitude + Latitude, data = ames_train)

wf = workflow() %>%
  add_recipe(rec) %>%
  add_model(krr_spec)

## 7.6 Hyperparameter tuning
cv_folds = vfold_cv(ames_train, v = 5, strata = Sale_Price)

param_grid = grid_regular(
  dials::penalty(range = c(-10, -3)), # lambda in {1e-10,...,1e-3}
  levels = 5
)

tune_results = tune_grid(
  wf,
  resamples = cv_folds,
  grid = param_grid,
  metrics = metric_set(rmse, rsq),
  control = control_grid(verbose = TRUE, save_pred = TRUE)
)
#autoplot(tune_results)


## 7.7 Model selection and final evaluation
best_params = select_best(tune_results, metric = "rmse")

final_spec = finalize_model(krr_spec, best_params)

final_wf = workflow() %>%
  add_recipe(rec) %>%
  add_model(final_spec)

final_fit   = final_wf %>% fit(data = ames_train)
pred = predict(final_fit, ames_test) %>%
  bind_cols(ames_test %>% select(Sale_Price))

metrics(pred, truth = Sale_Price, estimate = .pred)

end_time = Sys.time()
elapsed = end_time - start_time
elapsed