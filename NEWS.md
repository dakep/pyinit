# Version 1.1.1:
  * Fix numerical issue when using ATLAS BLAS libraries for CHOL decomposition.
    Manually check an approximation to the condition number to determine if the OLS
    coefficients can be computed.
  * Fix segfault on Oracle Solaris caused by an invalid call to lapack routine DSYEVR.

# Version 1.0.4:
  * Fix issues with uninitialized memory if OLS coefficients cannot be computed.

# Version 1.0.3:
  * Fix an issue with exact fits and if the M-scale estimate is (close to) 0.

# Version 1.0.2:
  * Initial release.
