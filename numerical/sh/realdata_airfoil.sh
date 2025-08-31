#!/bin/bash

mkdir -p output/log/realdata

for setting in cov_shift_null prior_shift_null cov_shift_alt prior_shift_alt; do
  echo "> Starting setting '${setting}'."

  name_output="output/log/realdata/airfoil__${setting}.output"
  Rscript numerical/realdata__airfoil.R $setting &> $name_output &
done

Rscript numerical/realdata__airfoil__local.R prior_shift_alt &> output/log/realdata/airfoil__local.output &

echo "> Waiting tasks to be finished."

wait

echo "> Finished. Please run 'extract_results.sh' to generate the tables."
