rm(list = ls())
start_time = Sys.time()
###### REQUIRED PACKAGES ----------------------------------------------
library(FastKRR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidymodels)
library(microbenchmark)

library(KRLS)
library(KRMM)
library(kernlab)


## 5. Structure and functionality of FastKRR package ----------------
## 5.2. Main functions for model fitting and prediction
#generate data ---------------------------------------------------------------
library(FastKRR)  
set.seed(1)
n = 1000 ; d = 1

# generating data
X = matrix(runif(n*d, 0, 1), nrow = n, ncol = d)
y = as.vector(sin(2*pi*X^3) + rnorm(n, 0, 0.1))


# grid for prediction
X_new = matrix(seq(0, 1, length.out = 1200), ncol = 1)
y_new = as.vector(sin(2*pi*X_new^3) + rnorm(1200, 0, 0.1))

#fit exact model ---------------------------------------------------------------
fit_exact = fastkrr(X, y, kernel = "gaussian", opt = "exact", verbose = FALSE)
pred_exact = predict(fit_exact, X_new)

FastKRR::error(fit_exact)
plot(fit_exact)



## 5-3. Interfaces for approximation methods
#approx examples ---------------------------------------------------------------
fit_nystrom = fastkrr(X, y, kernel = "gaussian",opt = "nystrom", m = 100, verbose = FALSE)
fit_pivoted = fastkrr(X, y, kernel = "gaussian",opt = "pivoted", m = 100, verbose = FALSE)
fit_rff = fastkrr(X, y, kernel = "gaussian",opt = "rff", m = 200, verbose = FALSE)

#krr-compare -------------------------------------------------------------------
set.seed(1)
n = 1000 ; d = 1 

# generating data
X = matrix(runif(n*d, 0, 1), nrow = n, ncol = d)
y = as.vector(sin(2*pi*X^3) + rnorm(n, 0, 0.1))

X_new = matrix(seq(0, 1, length.out = 1200), ncol = 1)
y_new = as.vector(sin(2*pi*X_new^3) + rnorm(1200, 0, 0.1))


fit_exact   = fastkrr(X, y, kernel = "gaussian", opt = "exact",   verbose = FALSE)
fit_nystrom = fastkrr(X, y, kernel = "gaussian", opt = "nystrom", m = 100, verbose = FALSE)
fit_pivoted = fastkrr(X, y, kernel = "gaussian", opt = "pivoted", m = 100, verbose = FALSE)
fit_rff     = fastkrr(X, y, kernel = "gaussian", opt = "rff",     m = 200, verbose = FALSE)


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




#approx_kernel_examples---------------------------------------------------------
K = make_kernel(X, kernel = "gaussian", rho = 1)
# approximation of kernel matrix
K_nystrom = approx_kernel(K = K, d = d, rho = 1, opt = "nystrom")
K_pchol = approx_kernel(K = K, d = d, rho = 1, opt = "pivoted")
K_rff = approx_kernel(X = X, d = d, rho = 1, kernel = "gaussian", opt = "rff")






## 5.5 Cross-validation and CVST support
#fit fastcv model --------------------------------------------------------------
exact_fast = fastkrr(X, y, kernel = "gaussian", opt = "exact", fastcv = TRUE)
nystrom_fast = fastkrr(X, y, kernel = "gaussian", opt = "nystrom", fastcv = TRUE)
pivoted_fast = fastkrr(X, y, kernel = "gaussian", opt = "pivoted", fastcv = TRUE)
rff_fast = fastkrr(X, y, kernel = "gaussian", opt = "rff", fastcv = TRUE)

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

pmse = function(model, exact) {
  X = attr(model, "x")
  n = nrow(X)
  
  n.new = ceiling(1.2 * n)
  rng = apply(X, 2, range)
  X.new = sapply(seq_len(ncol(X)), function(j) runif(n.new, rng[1, j], rng[2, j]))
  K.new = make_kernel(X, X.new, kernel = attr(model,"kernel"), rho = attr(model,"rho"))
  
  tab = c(0, 0)
  names(tab) = c("squared", "abs")
  
  if (attr(model, "opt") == "rff") {
    pred_model = predict(model, X.new)       
    pred_exact = K.new %*% attr(exact,"coefficients")
  } else {
    pred_model = K.new %*% attr(model,"coefficients")
    pred_exact = K.new %*% attr(exact,"coefficients")
  }
  
  tab["squared"] = mean((pred_exact - pred_model)^2)
  tab["abs"] = mean(abs(pred_exact - pred_model))
  return(tab)
}


pmse_plot = function(df, size= 12){
  ggplot(df, aes(x = expr, y = value)) +
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(~ metric, scales = "free_y") +
    labs(x = NULL, y = "PMSE") +
    theme_minimal() + theme(axis.text.x = element_text(size = size))
}





## 6.1 Computation time with external packages
#Computation-time-comparison ---------------------------------------------------
set.seed(1)
n = 1000
d = 3     
X = matrix(rnorm(n*d,0, 1),n,d)
y = as.vector(sin(2*pi*rowMeans(X)^3) + rnorm(n, 0, 0.1))
lambda  = 1e-4
rho = 1



library(kernlab)
krr = function(X, y, lambda, rho){
  n_row = nrow(X)
  K = kernelMatrix(rbfdot(sigma=0.5), X)
  alpha_hat = solve(K + lambda * diag(n_row)) %*% y
  pred_y = t(alpha_hat) %*% K
  return("pred_y" = pred_y)
}


library(KRLS)
library(KRMM)
krr_compare = microbenchmark::microbenchmark(
  KRLS = krls(X= X, y= as.vector(y), whichkernel = "gaussian",
              lambda = lambda, sigma = rho, derivative = FALSE),
  KRMM = Kernel_Ridge_MM(Y_train = y,
                         X_train = X, # random effect
                         Matrix_covariates_train = X, # making kernel
                         kernel = "Gaussian",
                         rate_decay_kernel = rho, # rho
                         method="RKHS"),
  kernlab = krr(X, y, lambda = lambda, rho = rho),
  FastKRR = fastkrr(X, y, kernel = "gaussian", opt = "exact", 
                    rho = rho, lambda = lambda, verbose = FALSE),
  times = 100,   
  unit = "ms"
)

save(krr_compare, file = "output/package_comp.Rdata")
df = as.data.frame(krr_compare)


summary_tbl = df %>%
  group_by(expr) %>%
  summarize(
    Mean   = mean(time)   / 1e9,
    SE     = sd(time) / sqrt(100) / 1e9,      
    Median = median(time) / 1e9,
    .groups = "drop"
  )

library(xtable)
xt = xtable(summary_tbl, digits = 3,align = c("c", rep("c", ncol(summary_tbl))))
print(xt, include.rownames = FALSE, floating = FALSE)





## 6.2 Computation time and approximation error in FastKRR
#simualtion---------------------------------------------------------------------
set.seed(1)
n = 5000 
d = 3
rho = 1 / d
lambda = 1e-10

X = matrix(runif(n*d, 0, 1), nrow = n, ncol = d)
y = as.vector(sin(2*pi*rowMeans(X)^3) + rnorm(n, 0, 0.1))
K = make_kernel(X)



compare_models = function(X, y, rho, lambda, B = 20, seed = 1) {
  if (!is.null(seed)) set.seed(seed)
  
  methods = c("exact", "nystrom", "pivoted", "rff")
  
  time_rows = list()
  pmse_rows = list()
  
  for (b in seq_len(B)) {
    models_iter = list()
    
    for (m in methods) {
      model_fit = NULL
      t_ns = microbenchmark(
        { model_fit = fastkrr(X, y,
                              kernel = "gaussian",
                              opt    = m,
                              rho    = rho,
                              lambda = lambda,
                              m = 200,  ##-------------------------------------------
                              eps = 1e-10, ##---------------------------------------
                              verbose = FALSE) },
        times = 1, unit = "ns"
      )$time
      
      time_rows[[length(time_rows) + 1]] = data.frame(
        expr = factor(m, levels = methods),
        time = as.numeric(t_ns)
      )
      
      models_iter[[m]] = model_fit
    }
    
    exact_model = models_iter[["exact"]]
    for (m in setdiff(methods, "exact")) {
      val = pmse(models_iter[[m]], exact_model)
      pmse_rows[[length(pmse_rows) + 1]] = data.frame(
        expr   = factor(m, levels = methods),
        metric = names(val),
        value  = as.numeric(val)
      )
    }
  }
  
  time_df = bind_rows(time_rows)
  pmse_df = bind_rows(pmse_rows)
  
  list(
    time_plot = time_plot(time_df),  
    pmse_plot = pmse_plot(pmse_df)   
  )
}



plot_n5000_d3_eps6 = compare_models(X, y, rho, lambda, B = 10)
save(plot_n5000_d3_eps6, file = "output/plot_n5000_d3_eps6.RData")


#time comparison ---------------------------------------------------------------
df_time = as.data.frame(plot_n5000_d3_eps6$time_plot$data)

summary_tbl = df_time %>%
  group_by(expr) %>%
  summarize(
    Mean   = mean(time)   / 1e9,
    SE = sd(time) /sqrt(100) /1e9,
    Median = median(time) / 1e9,
    .groups = "drop"
  )

colnames(summary_tbl)[colnames(summary_tbl) == "expr"] = "Option"
xt = xtable(summary_tbl, digits = 3,align = c("c", rep("c", ncol(summary_tbl))))
print(xt, include.rownames = FALSE, floating = FALSE)


#pmse comparison ---------------------------------------------------------------
df_pmse = as.data.frame(plot_n5000_d3_eps6$pmse_plot$data)

summary_pmse = df_pmse %>%
  mutate(Measure = ifelse(metric == "squared", "MSE", "MAE")) %>%
  group_by(Measure, expr) %>%
  summarize(
    Mean     = mean(value) * 1000,
    SE       = sd(value) / sqrt(100) * 1000,
    Median   = median(value) * 1000,
    .groups = "drop"
  )


summary_pmse$Measure = factor(summary_pmse$Measure, levels = c("MSE", "MAE"))
summary_pmse = summary_pmse %>% arrange(Measure, expr)
colnames(summary_pmse)[colnames(summary_pmse) == "expr"] = "Option"
summary_pmse = as.data.frame(summary_pmse)
summary_pmse$Measure = c("MSE", "", "", "MAE", "", "") 

xt = xtable(summary_pmse, digits = 5, align = c("c", rep("c", ncol(summary_pmse))))
addtorow = list(
  pos = list(3),  
  command = "\\hline \n"
)

print(xt,
      include.rownames = FALSE,
      floating = FALSE,
      add.to.row = addtorow)




## 6.3 Model selection via cross-validation 
#simualtion---------------------------------------------------------------------
set.seed(1)
n = 2000 
d = 3
rho = 1
lambda = seq(1e-5, 1e-12, length.out = 50)

X = matrix(runif(n*d, 0, 1), nrow = n, ncol = d)
y = as.vector(sin(2*pi*rowMeans(X)^3) + rnorm(n, 0, 0.1))
K = make_kernel(X)

exact_exact = fastkrr(X, y, kernel = "gaussian", opt = "exact", rho = rho, fastcv = FALSE, verbose = FALSE)
exact_fastCV = fastkrr(X, y, kernel = "gaussian", opt = "exact",lambda = lambda, rho = rho, fastcv = TRUE, verbose = FALSE)
benchmark_result_exact = microbenchmark::microbenchmark(
  "exact_exact" = fastkrr(X, y, kernel = "gaussian", opt = "exact", rho = rho, lambda = lambda,verbose = FALSE,fastcv = FALSE),
  "exact_fastCV" = fastkrr(X, y, kernel = "gaussian", opt = "exact", rho = rho, lambda = lambda,verbose = FALSE, fastcv = TRUE),
  times = 100 
)

save(exact_exact, file = "output/exact_exact.RData")
save(exact_fastCV, file = "output/exact_fastCV.RData")
save(benchmark_result_exact, file = "output/benchmark_result_exact.RData")
df = as.data.frame(benchmark_result_exact)


df$expr = factor(df$expr,
                 levels = c("exact_exact", "exact_fastCV"),
                 labels = c("Exact CV", "Fast CV"))

summary_cv = df %>%
  group_by(expr) %>%
  summarize(
    Mean     = mean(time) / 1e9,
    Variance = sd(time) /sqrt(100) / (1e9),
    Median   = median(time) / 1e9,
    .groups = "drop"
  )

xt = xtable(summary_cv, digits = 3, align = c("c", rep("c", ncol(summary_cv))))
print(xt, include.rownames = FALSE, floating = FALSE)






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

::metrics(pred, truth = Sale_Price, estimate = .pred)

sessionInfo()
end_time = Sys.time()
elapsed = end_time - start_time
elapsed