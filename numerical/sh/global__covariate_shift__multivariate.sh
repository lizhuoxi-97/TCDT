#!/bin/bash

num_simu=1000
p=4

setting_dgp=different_mean
dist_x=norm
dist_eps=t
shift_x=true

mkdir -p "output/log/simulations/global/covariate_shift/${setting_dgp}"

for setting_m in sindex_quad; do
  for n in 100 200; do
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

echo "> Waiting tasks to be finished."

wait

echo "> Finished. Please run 'extract_results.sh' to generate the tables and figures."