# TCDT: Two-Sample Conditional Distribution Testing

This R package, `TCDT`, provides tools for performing two-sample tests to determine 
if two conditional distributions are the same. 
It implements the methods proposed by Yan, Li and Zhang (2025), using either 
conditional energy distance (CED) or conditional maximum mean discrepancy (CMMD).

The package supports both global and local tests:
- **Global tests** assess the equality of conditional distributions across the 
entire support of the conditioning variables.
- **Local tests** focus on the equality of conditional distributions at a specific point.

## Installation

You can install the `TCDT` package from GitHub using the `devtools` package:

```R
# install.packages("devtools")
devtools::install_github("lizhuoxi-97/TCDT")
```

## Usage

The primary function in this package is `tcdt()`. 
Below are examples of how to use it for global and local testing.

### Global Test

A global test checks if the conditional distribution of `Y` given `X` is the same in two samples.

```R
# Load the package
library(TCDT)

# Generate some sample data
set.seed(42)
n1 <- 100
n2 <- 100
p <- 2
X1 <- matrix(rnorm(n1 * p), ncol = p)
X2 <- matrix(rnorm(n2 * p), ncol = p)
Y1 <- matrix(rnorm(n1), ncol = 1)
Y2 <- matrix(rnorm(n2), ncol = 1)

# Perform a global test
result_global <- tcdt(X1, X2, Y1, Y2, stat = "ced")

# Print the results
print(result_global)
```

### Local Test

A local test checks if the conditional distribution of `Y` given `X` is the same 
in two samples at a specific point `x0`.

```R
# Load the package
library(TCDT)

# Generate some sample data
set.seed(42)
n1 <- 100
n2 <- 100
p <- 2
X1 <- matrix(rnorm(n1 * p), ncol = p)
X2 <- matrix(rnorm(n2 * p), ncol = p)
Y1 <- matrix(rnorm(n1), ncol = 1)
Y2 <- matrix(rnorm(n2), ncol = 1)

# Specify the point for the local test
x0 <- c(0, 0)

# Perform a local test
result_local <- tcdt(X1, X2, Y1, Y2, x0 = x0, stat = "ced")

# Print the results
print(result_local)
```

## Reproducing Simulation Results

The `/numerical` directory contains the scripts to reproduce the simulation and 
realdata results from Yan, Li and Zhang (2025).

### Running Simulations and Realdata Examples

The shell scripts in `/numerical/sh` are used to run the simulations or realdata examples. 
Each script corresponds to a specific simulation example or realdata example described in the paper.

For example, to run the simulations for the local univariate case (Example 1 in the paper), 
you can execute the following command from the root of the project:

```bash
bash numerical/sh/local__univariate.sh
```

This will run the corresponding R script (`numerical/simulations__local__univariate.R`) 
with the appropriate settings, which might utilize multiple cores simultaneously.
One can explore other custom data settings of interest by modifying the specific 
R scripts under `numerical/`.

### Extracting Results

After running the simulations or realdata examples, the raw results are stored 
in the `/output/results` directory. 
To generate the tables and figures presented in the paper, you can run the 
`extract_results.sh` shell script located in the `/numerical/sh` directory.

```bash
bash numerical/sh/extract_results.sh
```

This script will process the raw results and generate the final tables and figures 
in the `/output/tables` and `/output/figures` directories, respectively.

## License

This project is licensed under the MIT License.

## References

Yan, J., Li, Z., \& Zhang, X. (2025). Distance and Kernel-Based Measures for 
Global and Local Two-Sample Conditional Distribution Testing. *arXiv preprint arXiv:2210.08149*.