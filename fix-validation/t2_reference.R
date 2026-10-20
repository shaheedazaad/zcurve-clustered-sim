# Independent check: maximise sum_i w_i log sum_k theta_k f_k(z_i) directly with optim()
source("fix-validation/common.R")
ref_fit <- function(fit_obj, w_x, w_b){
  ctl <- fit_obj$control; mu <- fit_obj$fit$mu
  z <- fit_obj$data; keep <- z > ctl$a & z < ctl$b
  cens <- fit_obj$data_censoring
  ll <- cbind(sapply(mu, function(m) zcurve:::.zdist_lpdf(z[keep], m, 1, ctl$a, ctl$b)))
  if(nrow(cens) > 0) ll <- rbind(ll, sapply(mu, function(m) mapply(function(l, u) zcurve:::.zdist_cens_lpdf(l, u, m, 1, ctl$a, ctl$b), cens$lb, cens$ub)))
  w <- c(w_x[keep], w_b)
  obj <- function(theta) sum(w * log(exp(ll) %*% theta))
  sm  <- function(par) { e <- exp(c(0, par)); e / sum(e) }
  best <- NULL
  for(r in 1:20){ o <- optim(rnorm(length(mu) - 1, 0, 2), function(p) -obj(sm(p)), method = "BFGS", control = list(maxit = 5000, reltol = 1e-14))
    if(is.null(best) || o$value < best$value) best <- o }
  prop_high <- (sum(w_x[z > ctl$b])) / (sum(w_x[z > ctl$a]) + sum(w_b))
  list(obj = obj, theta = sm(best$par), loglik = -best$value,
       est = zcurve:::.get_estimates(mu, sm(best$par), prop_high, ctl$sig_level, ctl$a))
}
weights_of <- function(f){ tab <- table(c(f$data_id, f$data_censoring$id)); list(x = 1 / as.vector(tab[as.character(f$data_id)]), b = 1 / as.vector(tab[as.character(f$data_censoring$id)])) }
report <- function(label, f){
  w <- weights_of(f); set.seed(99); r <- ref_fit(f, w$x, w$b)
  cat(sprintf("%-34s zcurve: ERR %.4f EDR %.4f logLik %.4f | reference optimum: ERR %.4f EDR %.4f logLik %.4f\n",
    label, f$coefficients["ERR"], f$coefficients["EDR"], r$obj(f$fit$weights), r$est["ERR"], r$est["EDR"], r$loglik))
}
# 1) the dataset from the GitHub issue (clusters of 1-5)
sim <- simulate_z_curve_dataset(target_err = 0.7, target_n_sig = 300, seed = simulation_seed, article_seed = article_seed)$data
report("issue dataset", fit_w(sim$z_string, sim$article_id))
# 2) the copied datasets from test 1, incl. censored + z>6 values
for(s in 1:3){
  sim <- simulate_z_curve_dataset(target_n_sig = 300, seed = 100 + s, mixture_info = mi, max_effects_per_article = 1)$data
  x <- c(sim$z_string, rep("p < 0.01", 30), "z = 6.5", "z = 7.2", "z = 8.0")
  set.seed(s); m <- sample(1:5, length(x), replace = TRUE)
  report(sprintf("copied 1-5x [seed %d]", s), fit_w(rep(x, m), rep(seq_along(x), m)))
}
