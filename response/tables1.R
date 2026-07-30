library(FastKRR)
library(KRMM)
library(dplyr)
library(xtable)
library(ggplot2)
library(microbenchmark)

set.seed(1)

time_plot = function(df, rm = 0){
  idx = df$expr %in% rm
  df = df[!idx, ]
  stats = df %>%
    group_by(expr) %>%
    summarize(Q1 = quantile(time/1e6, 0.25), Q3 = quantile(time/1e6, 0.75),
              IQR = Q3 - Q1, .groups = "drop")
  ymin = min(stats$Q1 - 1.5*stats$IQR); ymax = max(stats$Q3 + 1.5*stats$IQR)
  ggplot(df, aes(x = expr, y = time/1e6)) +
    geom_boxplot(outlier.shape = NA) +
    coord_cartesian(ylim = c(ymin, ymax)) +
    labs(x = "Method", y = "Time (milliseconds)") +
    theme_minimal()
}

calc_metric = function(pkg_error_list) {
  mse = sapply(pkg_error_list, function(e) mean(as.vector(e)^2, na.rm = TRUE))
  mae = sapply(pkg_error_list, function(e) mean(abs(as.vector(e)), na.rm = TRUE))
  data.frame(MSE_mean = mean(mse), MSE_se = sd(mse)/sqrt(length(mse)),
             MAE_mean = mean(mae), MAE_se = sd(mae)/sqrt(length(mae)))
}

make_tbl = function(time_df, pred_list, label) {
  t_tbl = time_df %>%
    group_by(expr) %>%
    summarize(Mean = mean(time)/1e9, SE = sd(time)/sqrt(n())/1e9,
              Median = median(time)/1e9, .groups = "drop") %>%
    mutate(Method = as.character(expr),
           `Mean (SE)` = sprintf("%.4f (%.4f)", Mean, SE),
           Median = sprintf("%.4f", Median)) %>%
    select(Method, `Mean (SE)`, Median)
  
  p_tbl = do.call(rbind, lapply(pred_list, calc_metric)) %>%
    tibble::rownames_to_column("Method") %>%
    mutate(`MSE (SE)` = sprintf("%.4f (%.4f)", MSE_mean, MSE_se),
           `MAE (SE)` = sprintf("%.4f (%.4f)", MAE_mean, MAE_se)) %>%
    select(Method, `MSE (SE)`, `MAE (SE)`)
  
  final = t_tbl %>% left_join(p_tbl, by = "Method")
  print(xtable(final, label = label, align = c("l","l","c","c","c","c")),
        include.rownames = FALSE, sanitize.text.function = identity,
        add.to.row = list(pos = list(-1),
                          command = paste0("\\hline\n",
                                           " & \\multicolumn{2}{c}{Time} & \\multicolumn{2}{c}{Prediction} \\\\\n",
                                           "\\cline{2-3} \\cline{4-5}\n")),
        hline.after = c(-1, 0, nrow(final)), caption.placement = "top")
  final
}

n = 1000
d = 3
rho = 1

gen_data = function(n, d) {
  X = matrix(runif(n * d, 0, 1), n, d)
  colnames(X) = paste0("X", seq_len(d))
  y = as.vector(sin(2 * pi * rowMeans(X)^3) + rnorm(n, 0, 0.1))
  list(X = X, y = y, data = data.frame(X, y = y))
}

train = gen_data(n, d)
X = train$X; y = train$y; data = train$data

n_z = ceiling(1.2 * n)
Z = matrix(runif(n_z * d, 0, 1), n_z, d)
colnames(Z) = paste0("X", seq_len(d))
y_z = as.vector(sin(2 * pi * rowMeans(Z)^3))
Z_df = as.data.frame(Z)

rate_decay_kernel = rho * d

exact_time = microbenchmark(
  KRMM = Kernel_Ridge_MM(Y_train = y, Matrix_covariates_train = X,
                         method = "RKHS", kernel = "Gaussian",
                         rate_decay_kernel = rate_decay_kernel),
  
  FastKRR_exact_REML = fastkrr(data = data, response = "y", kernel = "gaussian",
                               opt = "exact", rho = rho,
                               selection_method = "REML", verbose = FALSE),
  
  FastKRR_exact_fastCV = fastkrr(data = data, response = "y", kernel = "gaussian",
                                 opt = "exact", rho = rho,
                                 selection_method = "fastCV", verbose = FALSE),
  
  FastKRR_nystrom_fastCV = fastkrr(data = data, response = "y", kernel = "gaussian",
                                   opt = "nystrom", rho = rho,
                                   selection_method = "fastCV", verbose = FALSE),
  
  FastKRR_pivoted_fastCV = fastkrr(data = data, response = "y", kernel = "gaussian",
                                   opt = "pivoted", rho = rho,
                                   selection_method = "fastCV", verbose = FALSE),
  
  FastKRR_rff_fastCV = fastkrr(data = data, response = "y", kernel = "gaussian",
                               opt = "rff", rho = rho,
                               selection_method = "fastCV", verbose = FALSE),
  
  times = 100, unit = "ms"
)

dir.create("output", showWarnings = FALSE)
save(exact_time, file = "output/reml_vs_cvst_exact_time.Rdata")
time_plot(exact_time)

## ---- 예측 정확도 비교 (반복 시뮬레이션) ----
compare_exact_pred = list(KRMM = list(), FastKRR_exact_REML = list(),
                          FastKRR_exact_fastCV = list(), FastKRR_nystrom_fastCV = list(),
                          FastKRR_pivoted_fastCV = list(), FastKRR_rff_fastCV = list())

for (rep in 1:100) {
  tr = gen_data(n, d)
  Zr = matrix(runif(n_z * d, 0, 1), n_z, d); colnames(Zr) = paste0("X", seq_len(d))
  y_zr = as.vector(sin(2 * pi * rowMeans(Zr)^3))
  
  m_krmm = Kernel_Ridge_MM(Y_train = tr$y, Matrix_covariates_train = tr$X,
                           method = "RKHS", kernel = "Gaussian",
                           rate_decay_kernel = rate_decay_kernel)
  pred_krmm = as.vector(Predict_kernel_Ridge_MM(m_krmm, Matrix_covariates_target = Zr))
  
  m_exact_reml    = fastkrr(data = tr$data, response = "y", kernel = "gaussian",
                            opt = "exact", rho = rho, selection_method = "REML",    verbose = FALSE)
  m_exact_fastcv  = fastkrr(data = tr$data, response = "y", kernel = "gaussian",
                            opt = "exact", rho = rho, selection_method = "fastCV",  verbose = FALSE)
  m_nystrom_fastcv = fastkrr(data = tr$data, response = "y", kernel = "gaussian",
                             opt = "nystrom", rho = rho, selection_method = "fastCV", verbose = FALSE)
  m_pivoted_fastcv = fastkrr(data = tr$data, response = "y", kernel = "gaussian",
                             opt = "pivoted", rho = rho, selection_method = "fastCV", verbose = FALSE)
  m_rff_fastcv     = fastkrr(data = tr$data, response = "y", kernel = "gaussian",
                             opt = "rff", rho = rho, selection_method = "fastCV", verbose = FALSE)
  
  compare_exact_pred$KRMM[[rep]]                    = y_zr - pred_krmm
  compare_exact_pred$FastKRR_exact_REML[[rep]]       = y_zr - predict(m_exact_reml,     newdata = as.matrix(Zr))
  compare_exact_pred$FastKRR_exact_fastCV[[rep]]     = y_zr - predict(m_exact_fastcv,   newdata = as.matrix(Zr))
  compare_exact_pred$FastKRR_nystrom_fastCV[[rep]]   = y_zr - predict(m_nystrom_fastcv, newdata = as.matrix(Zr))
  compare_exact_pred$FastKRR_pivoted_fastCV[[rep]]   = y_zr - predict(m_pivoted_fastcv, newdata = as.matrix(Zr))
  compare_exact_pred$FastKRR_rff_fastCV[[rep]]       = y_zr - predict(m_rff_fastcv,     newdata = as.matrix(Zr))
}
save(compare_exact_pred, file = "output/reml_vs_cvst_exact_pred.Rdata")

tbl_exact = make_tbl(exact_time, compare_exact_pred, "tab:reml-vs-cvst-exact")
tbl_exact


