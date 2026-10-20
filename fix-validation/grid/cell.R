# One simulation cell. Usage: Rscript cell.R <lib> <ci|pt> <dep> <clus> <err> <n> <reps>
#   dep : none | shared | hetero | corr       clus: 1-3 | 1-5 | 1-10 | skew       err: target ERR of the 0/2/4 mixture
# Each article draws a true effect, runs k tests, and reports only its significant ones.
args <- commandArgs(TRUE)
lib <- args[1]; do_ci <- args[2] == "ci"; dep <- args[3]; clus <- args[4]
err <- as.numeric(args[5]); n_target <- as.integer(args[6]); reps <- as.integer(args[7])
suppressMessages({library(zcurve, lib.loc = lib); library(tibble); library(dplyr)})
src <- readLines("fix-validation/sim_functions.R"); eval(parse(text = src[!grepl("^library", src)]))
mi  <- find_mixture_weights(err); ncp <- c(0, 2, 4); a <- qnorm(.975)
g   <- function(mu){ p1 <- 1 - pnorm(a, mu); pn <- pnorm(-a, mu); (p1^2 + pn^2) / (p1 + pn) }  # ERR of one result
rk  <- function() switch(clus, "1-3" = sample(1:3, 1), "1-5" = sample(1:5, 1), "1-10" = sample(1:10, 1),
                         "skew" = min(20, 1 + rgeom(1, 0.3)))
gen <- function(seed){
  set.seed(seed); rows <- list(); n <- 0; j <- 0
  while(n < n_target){
    j <- j + 1; k <- rk(); mu_i <- sample(ncp, 1, prob = mi$mixture_weights)
    mu <- switch(dep, none = sample(ncp, k, TRUE, mi$mixture_weights), shared = rep(mu_i, k),
                 hetero = mu_i + rnorm(k, 0, 0.5), corr = rep(mu_i, k))
    z  <- if(dep == "corr") mu + sqrt(.8) * rnorm(1) + sqrt(.2) * rnorm(k) else rnorm(k, mu)
    sig <- abs(z) > a
    if(any(sig)){ rows[[j]] <- data.frame(id = j, z = abs(z[sig]), mu = mu[sig]); n <- n + sum(sig) }
  }
  d <- do.call(rbind, rows)
  if(dep == "none") d$id <- sample(d$id)   # ids unrelated to effects
  d
}
one <- function(r){
  seed <- 1e6 * match(dep, c("none","shared","hetero","corr")) + 1e5 * match(clus, c("1-3","1-5","1-10","skew")) + 1e4 * round(err * 10) + n_target + r * 7919
  d <- gen(seed); x <- sprintf("z=%.8f", d$z); dat <- zcurve_data(x, id = d$id, rounded = FALSE)
  out <- c(rep = r, n = nrow(d), n_art = length(unique(d$id)), true_res = mean(g(d$mu)), true_art = mean(tapply(g(d$mu), d$id, mean)))
  res <- tryCatch({
    set.seed(r); w <- suppressWarnings(zcurve_clustered(dat, method = "w", bootstrap = if(do_ci) 200 else FALSE))
    o <- c(w = unname(w$coefficients["ERR"]))
    if(do_ci){
      set.seed(r); u <- zcurve(data = zcurve_data(x, rounded = FALSE), bootstrap = FALSE)
      set.seed(r); b <- suppressWarnings(zcurve_clustered(dat, method = "b", bootstrap = 100))
      cw <- summary(w)$coefficients["ERR", ]; cb <- summary(b)$coefficients["ERR", ]
      o <- c(o, w_l = cw[["l.CI"]], w_u = cw[["u.CI"]], u = unname(u$coefficients["ERR"]),
             b = cb[["Estimate"]], b_l = cb[["l.CI"]], b_u = cb[["u.CI"]])
    }
    o
  }, error = function(e) c(w = NA))
  c(out, res)
}
res <- do.call(dplyr::bind_rows, lapply(parallel::mclapply(seq_len(reps), one, mc.cores = 7), function(v) as.data.frame(t(v))))
res <- cbind(dep = dep, clus = clus, err = err, n_target = n_target, res)
write.csv(res, sprintf("fix-validation/grid/out/%s_%s_%s_%s_%d_%s.csv", basename(lib), dep, clus, err, n_target, args[2]), row.names = FALSE)
