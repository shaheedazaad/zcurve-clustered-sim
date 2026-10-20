# Validation of the fix for zcurve issue #33

Scripts that compare the original weighted clustered z-curve (`zcurve_clustered(method = "w")`) with the fixed version proposed in [FBartos/zcurve#33](https://github.com/FBartos/zcurve/issues/33).

## Setup

Install the two builds of zcurve into separate library folders, named `original` and `fixed` (the result files are named after them):

```sh
git clone https://github.com/FBartos/zcurve.git
mkdir -p libs/original libs/fixed
R CMD INSTALL --library=libs/original zcurve                       # current master
git -C zcurve fetch https://github.com/shaheedazaad/zcurve.git fix-weighted-clustered-em
git -C zcurve checkout FETCH_HEAD
R CMD INSTALL --library=libs/fixed zcurve                          # with the fix
```

The scripts also need `tibble` and `dplyr`. Run every script from the root of this repository, passing the library folder as the first argument.

## Scripts

| Script | What it checks |
|---|---|
| `t1_invariance.R <lib>` | Cases with a known answer, including censored p-values and z > 6: (1) all clusters of equal size, where weighting must not change the fit; (2) every result copied 1–5 times inside its own cluster with weight 1/copies, where the fit must match the uncopied data. |
| `t2_reference.R <lib>` | Maximises the weighted log-likelihood `sum_i w_i log sum_k theta_k f_k(z_i)` directly with `optim()` and compares it with the zcurve fit. |
| `grid/run_all.sh` | The simulation grid: 4 dependence structures × 4 cluster-size distributions × 3 true ERRs at 300 significant results, plus sample sizes 100 and 1000. 200 datasets per cell. Takes about 4 hours on 7 cores. Set `LIB_ORIG`, `LIB_FIXED` and `LIB_DEPS` first. |
| `grid/summarise.R` | Bias, RMSE, and 95% CI coverage and width per cell, written to `grid/summary.csv`. |

`sim_functions.R` contains the data-generating functions from `example.qmd`, based on Ulrich Schimmack's z-curve 3.0 tutorial. `common.R` loads the chosen zcurve build and holds the fitting helpers.

## Simulation design

Each article draws a true effect from a mixture of non-centrality values 0, 2 and 4 (weights chosen to hit the target ERR). It runs k tests and reports only the significant ones.

- **Dependence:**
  - `none`: independent results with random article IDs.
  - `shared`: all tests in an article share one effect.
  - `hetero`: tests vary around the article's effect with SD 0.5.
  - `corr`: one shared effect plus correlated noise, ρ = 0.8, which is close to near-duplicate analyses.
- **Tests per article (k):** uniform on 1–3, 1–5 or 1–10, or `skew` = 1 + geometric(0.3), capped at 20.

Estimates are compared with the **article-level ERR**, where each article counts once (the target of both clustering adjustments). Results also include the result-level ERR, where each reported result counts once.
