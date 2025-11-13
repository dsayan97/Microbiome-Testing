Microbiome-Testing
==================

This repository contains R code accompanying the article:

Bhattacharjee et al. (2024).
“Testing Microbiome Community Differences in High Dimensions:
A Bootstrap Approach for Compositional Data.”

The code implements:

- the proposed bootstrap two-sample test (EBC.R),
- the MANOVA extension (KEBC.R),
- the benchmark test of Cao, Lin & Li (2017) (CLL.R), and
- the simulation scripts used in Section 6 of the paper.

All functions are written in base R. Run_TwoSample.R relies only on the `highmean` package for data-generation utilities, and Table_1.R and Table_4.R rely on `parallel`, `foreach`, and `doParallel` for parallel computation.

------------------------------------------------------
1. Repository Structure
------------------------------------------------------

EBC.R               – Proposed two-sample bootstrap test

CLL.R               – Cao, Lin & Li (2017) test

KEBC.R              – Proposed MANOVA (K-sample)  test

Parameters_Two.R    – Parameters for two-sample simulations

Parameters_MANOVA.R – Parameters for MANOVA simulations

Run_TwoSample.R     – Script to produce values in Tables 1, 2, and 3 in the main manuscript (values might differ due to randomness)

Run_MANOVA.R        – Script to produce values in Table 4 in the main manuscript (values might differ due to randomness)

Table_1.R           - Script to reproduce Table 1 in the main manuscript (use `set.seed(1)` to match the values in the main manuscript)

Table_4.R           - Script to reproduce Table 4 in the main manuscript (use `set.seed(1)` to match the values in the main manuscript)

------------------------------------------------------
2. Core Functions
------------------------------------------------------

2.1 EBC.R — Proposed Two-Sample Test
------------------------------------

Function:
    EBC(X, Y, N = 1000)

Inputs:
- X, Y: Data matrices of size $n_1 \times p$ and $n_2 \times p$ (sample size x dimension)
- N: Number of bootstrap replicates (default: 1000)

Output:
- A bootstrap p-value for testing equality of mean vectors.

------------------------------------

2.2 CLL.R — Cao, Lin & Li (2017) Test
------------------------------------

Function:
    CLL(X, Y, alpha = 0.05)

Inputs:
- X, Y: Data matrices ($n \times p$) after CLR transformation (see equation 3 in the main manuscript)
- alpha: Significance level (default: 0.5)

Output:
- 1 = reject H0
- 0 = fail to reject H0

------------------------------------

2.3 KEBC.R — Proposed MANOVA-Type Test (K-Sample)
--------------------------------------------------

Function:
    KEBC(X, N = 1000)

Inputs:
- X: A list of data matrices from K populations after CLR transformation (see equation 3 in the main manuscript)

    Example: X = list(Y1, Y2, ..., YK)

    Each Yk is $n_k \times p$, and sample sizes ($n_k$) may differ.
- N: Number of bootstrap samples.

Output:
- A bootstrap p-value for testing equality of mean vectors across K groups.

------------------------------------------------------
3. Parameter Files
------------------------------------------------------

Parameters_Two.R
    contains data generation codes for the two-sample simulations.

Parameters_MANOVA.R
    contains analogous settings for multigroup simulations.

------------------------------------------------------
4. Running the Simulation Scripts
------------------------------------------------------

4.1 Two-Sample Simulations — Run_TwoSample.R
---------------------------------------------

This script reproduces the results in Tables 1, 2 and 3.

Before running install library
    `highmean` and keep
    "EBC.R",
    "CLL.R",
    and "Parameters_Two.R" in the working directory.

Important parameters:
- n1, n2: Sample sizes for the two groups
- p: Dimension
- s: Sparsity level of the mean vector
- d: Separation parameter ($\delta$)
- alpha: Significance level
- ITER: Number of Monte Carlo iterations
- COV: Covariance structure (1 = (M.1), 2 = (M.2), 3 = (M.3))
- SETUP: Distribution (1 = Gaussian, 2 = t_5)

The example settings produce fourth, fifth and sixth values in row 10 of Table 1 for EBC, CLL, and MECAF.

---------------------------------------------

4.2 Multigroup Simulations — Run_MANOVA.R
---------------------------------------------

Reproduces the results in Table 4.

Before running keep
    "KEBC.R",
    and "Parameters_MANOVA.R" in the working directory.

Important parameters:
- n1–n5: Sample sizes for up to five populations
- p, s, d, alpha, ITER, SETUP: Same meaning as above

The example settings produce fourth, fifth and sixth values in row 10 of Table 4 for for K = 3, 4, and 5 populations.

------------------------------------------------------
5. Quick Start Examples
------------------------------------------------------

Two-Sample Test:
----------------
    G <- diag(10) - matrix(1, 10, 10)/10
    X <- matrix(rnorm(300), nrow = 30); X <- X %*% G
    Y <- matrix(rnorm(200), nrow = 20); Y <- Y %*% G

    pval <- EBC(X, Y, N = 500)
    print(pval)

K-Sample MANOVA Test:
---------------------
    G <- diag(10) - matrix(1, 10, 10)/10
    Y1 <- matrix(rnorm(300), nrow = 30); Y1 <- Y1 %*% G
    Y2 <- matrix(rnorm(250), nrow = 25); Y2 <- Y2 %*% G
    Y3 <- matrix(rnorm(200), nrow = 20); Y3 <- Y3 %*% G

    pval <- KEBC(list(Y1, Y2, Y3), N = 500)
    print(pval)

------------------------------------------------------
6. Notes on Reproducibility
------------------------------------------------------

- Run Table_1.R and Table_4.R with `set.seed(1)` to reproduce tables 1 and 4 in the main manuscript.
- Run Table_1.R with `cov.sim = 2`(line 13) and `cov.sim = 3`(line 13) with `set.seed(1)` to reproduce tables 2 and 3, respectively, in the main manuscript.
- Warning!! Compilation time for Table_1.R and Table_4.R can be really long. Instead, use Parameters_Two.R and Parameters_MANOVA.R to produce parts of the tables 1-4. However, in this case, values might *not exactly match* due to randomness.
- All data matrices should be CLR transformed (rows sum to 0).
- Bootstrap-based methods may be computationally intensive; reduce N for testing if necessary.

------------------------------------------------------
7. Contact
------------------------------------------------------

For questions or assistance, please contact the corresponding author listed in the paper.
