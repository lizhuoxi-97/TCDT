#!/bin/bash

setting_dgp="different_mean"
num_simu=1000
shift_y=true

mkdir -p "output/log/simulations/global/label_shift/${setting_dgp}"

# Need 3 * 2 * 1 * 6 = 36 kernels to run simultaneously
for n in 50 100 200; do
  for p in 5 20; do
    for setting_m in nonlinear; do
      for signal_m in 0 0.2 0.4 0.6 0.8 1; do
        name_output="output/log/simulations/global/label_shift/${setting_dgp}/"\
"m_${setting_m}__"\
"${signal_m}__"\
"n_${n}_${n}__"\
"p_${p}__"\
"shift_y_${shift_y}"\
".output"

        Rscript numerical/simulations__global__label_shift.R \
        $setting_dgp $n $n $p $shift_y $setting_m $signal_m $num_simu &> $name_output &
      done
    done
  done
done

echo "> Waiting tasks to be finished."

wait

echo "> Finished. Please run 'extract_results.sh' to generate the tables and figures."