#!/bin/bash

mkdir -p "output/log/simulations/local/multivariate"
name_output="output/log/simulations/local/multivariate/log.output"
Rscript numerical/simulations__local__multivariate.R &> $name_output &

echo "> Waiting tasks to be finished."

wait

echo "> Finished. Figures have also been generated."
