# Summarise every finished cell of the grid into one table (fix-validation/grid/summary.csv)
suppressMessages(library(dplyr))
fx <- list.files("fix-validation/grid/out", "^fixed_.*_ci\\.csv$", full.names = TRUE)
og <- list.files("fix-validation/grid/out", "^original_.*_pt\\.csv$", full.names = TRUE)
f <- bind_rows(lapply(fx, read.csv, colClasses = c(clus = "character")))
o <- bind_rows(lapply(og, read.csv, colClasses = c(clus = "character"))) |> select(dep, clus, err, n_target, rep, w_orig = w)
d <- left_join(f, o, by = c("dep", "clus", "err", "n_target", "rep"))
cov <- function(l, u, t) mean(l <= t & t <= u, na.rm = TRUE)
s <- d |> group_by(dep, clus, err, n_target) |> summarise(
  reps = n(), fails = sum(is.na(w)), n_art = mean(n_art),
  true_res = mean(true_res), true_art = mean(true_art),
  bias_u = mean(u - true_art, na.rm = TRUE), bias_w_orig = mean(w_orig - true_art, na.rm = TRUE),
  bias_w = mean(w - true_art, na.rm = TRUE), bias_b = mean(b - true_art, na.rm = TRUE),
  rmse_u = sqrt(mean((u - true_art)^2, na.rm = TRUE)), rmse_w_orig = sqrt(mean((w_orig - true_art)^2, na.rm = TRUE)),
  rmse_w = sqrt(mean((w - true_art)^2, na.rm = TRUE)), rmse_b = sqrt(mean((b - true_art)^2, na.rm = TRUE)),
  cov_w = cov(w_l, w_u, true_art), cov_b = cov(b_l, b_u, true_art),
  width_w = mean(w_u - w_l, na.rm = TRUE), width_b = mean(b_u - b_l, na.rm = TRUE), .groups = "drop") |>
  mutate(dep = factor(dep, c("none", "shared", "hetero", "corr")), clus = factor(clus, c("1-3", "1-5", "1-10", "skew"))) |>
  arrange(dep, n_target != 300, clus, err, n_target)
write.csv(s, "fix-validation/grid/summary.csv", row.names = FALSE)
options(width = 250)
print(as.data.frame(s |> mutate(across(where(is.double), ~ round(.x, 3)))), row.names = FALSE)
