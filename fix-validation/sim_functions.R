library(tibble); library(dplyr)
library(zcurve)

target_n <- 300
target_err <- 0.70
simulation_seed <- 20260521
article_seed <- 20260522
cluster_bootstrap <- 500
parallel <- TRUE

article_min_effects <- 1
article_max_effects <- 5

# From Ulrich Schimmack's z-curve 3.0 tutorial

find_mixture_weights <- function(target_err = 0.70,
                                 noncentrality_values = c(0, 2, 4),
                                 weight_step = 0.01) {
  z_critical <- qnorm(0.975)
  directional_power <- pnorm(abs(noncentrality_values), mean = z_critical)
  sign_error_rate <- 1 - pnorm(noncentrality_values, mean = -z_critical)
  significance_rate <- directional_power + sign_error_rate
  extreme_tail_rate <- pnorm(6, mean = noncentrality_values, lower.tail = FALSE)

  candidate_weights <- expand.grid(rep(list(seq(0, 1, by = weight_step)),
                                       length(noncentrality_values)))
  names(candidate_weights) <- paste0("weight_", noncentrality_values)

  candidate_weights <- candidate_weights[
    abs(rowSums(candidate_weights) - 1) < 1e-12,
  ]

  candidate_weights$err <- apply(
    candidate_weights,
    1,
    \(mixture_weights) {
      extreme_weight <- sum(mixture_weights * extreme_tail_rate)
      observed_weights <- c(
        mixture_weights * (1 - extreme_tail_rate),
        extreme_weight
      )
      observed_significance <- c(significance_rate, 1)
      significant_weights <- observed_weights * observed_significance
      significant_weights <- significant_weights / sum(significant_weights)

      sum(significant_weights * c(directional_power, 1))
    }
  )
  candidate_weights$distance_to_target <- abs(candidate_weights$err - target_err)

  best_weights <- candidate_weights[which.min(candidate_weights$distance_to_target), ]
  mixture_weights <- as.numeric(best_weights[seq_along(noncentrality_values)])

  list(
    mixture_weights = setNames(
      mixture_weights,
      paste0("ncp_", noncentrality_values)
    ),
    approximate_err = best_weights$err,
    expected_sig_rate = sum(mixture_weights * significance_rate)
  )
}


simulate_z_curve_dataset <- function(target_err = 0.70,
                                     target_n_sig = 500,
                                     noncentrality_values = c(0, 2, 4),
                                     weight_step = 0.01,
                                     mixture_info = NULL,
                                     seed = 123,
                                     min_effects_per_article = 1,
                                     max_effects_per_article = 5,
                                     article_seed = NULL) {
  set.seed(seed)

  z_critical <- qnorm(0.975)
  if (is.null(mixture_info)) {
    mixture_info <- find_mixture_weights(
      target_err = target_err,
      noncentrality_values = noncentrality_values,
      weight_step = weight_step
    )
  }

  mixture_weights <- unname(mixture_info$mixture_weights)
  expected_sig_rate <- mixture_info$expected_sig_rate

  significant_z_scores <- numeric(0)
  total_draws <- 0L

  while (length(significant_z_scores) < target_n_sig) {
    remaining_n_sig <- target_n_sig - length(significant_z_scores)
    batch_size <- max(
      500L,
      ceiling(remaining_n_sig / expected_sig_rate * 1.25)
    )

    component_id <- sample(
      x = seq_along(noncentrality_values),
      size = batch_size,
      replace = TRUE,
      prob = mixture_weights
    )

    z_scores <- rnorm(
      n = batch_size,
      mean = noncentrality_values[component_id],
      sd = 1
    )

    significant_z_scores <- c(
      significant_z_scores,
      abs(z_scores[abs(z_scores) > z_critical])
    )
    total_draws <- total_draws + batch_size
  }

  significant_z_scores <- significant_z_scores[seq_len(target_n_sig)]

  if (!is.null(article_seed)) {
    set.seed(article_seed)
  }

  article_sizes <- integer(0)
  effects_remaining <- target_n_sig

  while (effects_remaining > 0) {
    next_article_size <- sample(
      min_effects_per_article:min(max_effects_per_article, effects_remaining),
      size = 1
    )

    article_sizes <- c(article_sizes, next_article_size)
    effects_remaining <- effects_remaining - next_article_size
  }

  article_ids <- sample(
    sprintf(
      "article_%03d",
      rep(seq_along(article_sizes), times = article_sizes)
    ),
    size = target_n_sig
  )

  simulated_dataset <- tibble(
    z = significant_z_scores,
    target_err = target_err,
    approximate_err = mixture_info$approximate_err,
    article_id = article_ids,
    z_string = sprintf("z=%.8f", z)
  )

  invisible(list(
    data = simulated_dataset,
    target_err = target_err,
    target_n_sig = target_n_sig,
    total_draws = total_draws,
    approximate_err = mixture_info$approximate_err,
    mixture_weights = setNames(
      mixture_weights,
      paste0("ncp_", noncentrality_values)
    )
  ))
}

# One Dataset With 1-5 Effects Per Article

