for example in local__univariate local__multivariate global__covariate_shift__univariate global__covariate_shift__multivariate global__label_shift realdata__airfoil; do
  echo "> Extracting results for ${example}..."
  Rscript numerical/extract_results.R $example
done
