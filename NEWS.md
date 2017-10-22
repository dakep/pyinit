# Version 0.2.5:
    * Ignore errors from eigenvalue decomposition related to convergence issues.

# Version 0.2.3:
    * Add better handling of singular design matrices. If the design matrix X is singular,
      the inverse of the square root of X is computed using the Moore-Penrose pseudo inverse.
      This affects the OLS estimates as well as the PSCs.
    * Fix bug in handling the intercept in design matrices with single covariates.

# Version 0.2.1:
    * Change handling of singular subsets. If a subset (obtained by removing observations based
      on a PSC and ordering) is singular, this subset is skipped instead of raising an error.
