#!/bin/bash

setting_dgp=different_mean
n=150
p=1
shift_x=true
dist_x=normt
setting=x
mode_local=parallel
num_simu=1000

mkdir -p "output/log/simulations/local/univariate/${setting_dgp}"

# Need 27 kernels to run simultaneously
for x0 in $(seq -1.0 0.1 1.0); do
  name_output="output/log/simulations/local/univariate/${setting_dgp}/"\
"n_${n}_${n}__"\
"p_${p}__"\
"setting_${setting}__"\
"dist_x_${dist_x}__"\
"shift_x_${shift_x}__"\
"mode_local_${mode_local}__"\
"x0_${x0}__"\
"num_simu_${num_simu}"\
".output"

  Rscript numerical/simulations__local__univariate.R \
  $setting_dgp $n $n $p $setting $dist_x $shift_x $mode_local $x0 $num_simu &> $name_output &
done

echo "> Waiting tasks to be finished."

wait

echo "> Finished. Please run 'extract_results.sh' to generate the tables and figures."
