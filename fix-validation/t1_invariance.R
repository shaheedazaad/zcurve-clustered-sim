# Exact checks: answers that a correct weighted likelihood must reproduce
source("fix-validation/common.R")
cmp <- function(label, a, b) cat(sprintf("%-48s ERR %.4f vs %.4f | EDR %.4f vs %.4f | max|dTheta| %.4f\n", label,
  a$coefficients["ERR"], b$coefficients["ERR"], a$coefficients["EDR"], b$coefficients["EDR"], max(abs(a$fit$weights - b$fit$weights))))
for(s in 1:3){
  sim <- simulate_z_curve_dataset(target_n_sig = 300, seed = 100 + s, mixture_info = mi, max_effects_per_article = 1)$data
  # add censored p-values and z > 6 values so every weighted code path is used
  x <- c(sim$z_string, rep("p < 0.01", 30), "z = 6.5", "z = 7.2", "z = 8.0")
  u <- fit_u(x)
  set.seed(s); id_eq <- sample(rep(seq_len(length(x) / 3), each = 3))
  cmp(sprintf("[seed %d] all clusters size 3  vs unclustered", s), fit_w(x, id_eq), u)
  set.seed(s); m <- sample(1:5, length(x), replace = TRUE)
  cmp(sprintf("[seed %d] copied 1-5x, w=1/m   vs unclustered", s), fit_w(rep(x, m), rep(seq_along(x), m)), u)
}
