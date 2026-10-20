#!/bin/zsh
# Runs the full grid: fixed build with CIs, original build point estimates only.
# Requires LIB_ORIG, LIB_FIXED (library dirs with each zcurve build) and LIB_DEPS (dependencies).
cd "$(dirname "$0")/../.."  # repo root
cells=()
for dep in none shared hetero corr; do
  for clus in 1-3 1-5 1-10 skew; do for err in 0.4 0.7 0.85; do cells+=("$dep $clus $err 300"); done; done
  for n in 100 1000; do cells+=("$dep 1-5 0.7 $n"); done
done
i=0
for c in "${cells[@]}"; do
  i=$((i+1)); set -- ${=c}
  f=fix-validation/grid/out/$(basename $LIB_FIXED)_$1_$2_$3_$4_ci.csv
  [[ -f $f ]] || R_LIBS=$LIB_DEPS Rscript fix-validation/grid/cell.R $LIB_FIXED ci $1 $2 $3 $4 200 >/dev/null 2>>fix-validation/grid/errors.log
  f=fix-validation/grid/out/$(basename $LIB_ORIG)_$1_$2_$3_$4_pt.csv
  [[ -f $f ]] || R_LIBS=$LIB_DEPS Rscript fix-validation/grid/cell.R $LIB_ORIG pt $1 $2 $3 $4 200 >/dev/null 2>>fix-validation/grid/errors.log
  echo "$(date +%H:%M:%S) done $i/${#cells[@]}: $c"
done
echo ALL DONE
