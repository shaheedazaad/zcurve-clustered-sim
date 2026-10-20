lib <- commandArgs(TRUE)[1]
suppressMessages({library(zcurve, lib.loc = lib); library(tibble); library(dplyr)})
src <- readLines("fix-validation/sim_functions.R"); eval(parse(text = src[!grepl("^library", src)]))
mi  <- find_mixture_weights(0.7)
fit_u <- function(strings, seed = 1){ set.seed(seed); zcurve(data = zcurve_data(strings, rounded = FALSE), bootstrap = FALSE) }
fit_w <- function(strings, id, seed = 1){ set.seed(seed); suppressWarnings(zcurve_clustered(zcurve_data(strings, id = id, rounded = FALSE), method = "w", bootstrap = FALSE)) }
