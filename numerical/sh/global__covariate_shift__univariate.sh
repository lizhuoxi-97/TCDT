#!/bin/bash

num_simu=1000
p=1

# Difference in conditional means -------------------------------

setting_dgp="different_mean"
dist_x=normt
shift_x=true
dist_eps=norm

mkdir -p "output/log/simulations/global/covariate_shift/${setting_dgp}"

# Need 2 * 2 * 6 = 24 kernels to run simultaneously
for n in 50 100; do
  for setting_m in "x2+x(x+2)(x-2)" "exp+sin"; do
    for signal_m in 0 0.4 0.8 1.2 1.6 2; do
      name_output="output/log/simulations/global/covariate_shift/${setting_dgp}/"\
"m_${setting_m}__"\
"${signal_m}__"\
"n_${n}_${n}__"\
"p_${p}__"\
"dist_x_${dist_x}__"\
"shift_x_${shift_x}__"\
"dist_eps_${dist_eps}"\
".output"

      Rscript numerical/simulations__global__covariate_shift.R \
      $setting_dgp $n $n $p $dist_x $shift_x $dist_eps $setting_m $signal_m $num_simu \
      &> $name_output &
    done
  done
done

echo "> Waiting tasks (diffm) to be finished."

# Difference in conditional variances -------------------------------

setting_dgp="different_var"
dist_x=normt
shift_x=true
dist_eps=norm
setting_m=quad

mkdir -p "output/log/simulations/global/covariate_shift/${setting_dgp}"

# Need 2 * 2 * 6 = 24 kernels to run simultaneously
for n in 100 200; do
  for setting_v in homo hetero; do
    for signal_v in 0 0.2 0.4 0.6 0.8 1; do
      name_output="output/log/simulations/global/covariate_shift/${setting_dgp}/"\
"v_${setting_v}__"\
"${signal_v}__"\
"m_${setting_m}__"\
"n_${n}_${n}__"\
"p_${p}__"\
"dist_x_${dist_x}__"\
"shift_x_${shift_x}__"\
"dist_eps_${dist_eps}"\
".output"

      Rscript numerical/simulations__global__covariate_shift.R \
      $setting_dgp $n $n $p $dist_x $shift_x $dist_eps $setting_m $setting_v $signal_v $num_simu \
      &> $name_output &
    done
  done
done

echo "> Waiting tasks (diffv) to be finished."

wait

echo "> Finished. Please run 'extract_results.sh' to generate the tables and figures."
