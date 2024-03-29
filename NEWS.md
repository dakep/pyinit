# Version 1.1.3:
  * Update BLAS includes to fix errors with calls to Fortran routines.

# Version 1.1.2:
  * Fix C code to import `<R_ext/Error.h>` manually. 
    Thanks to Prof. Ripley for providing the fix alongside the error report.
  * Address implicit conversion warnings on M1 macs by making the conversions explicit.

# Version 1.1.1:
  * Estimates from ill-conditioned subsets of the data are not returned anymore.
    Any estimate with a value of the objective function of 0 is due to an exact fit.
  * Fix numerical issue when using ATLAS BLAS libraries for CHOL decomposition.
    pyinit now manually checks if the reciprocal condition number is large enough to
    compute the OLS coefficients. If not, the subset is considered ill-conditioned.
  * Fix segfault on Oracle Solaris caused by an invalid call to lapack routine DSYEVR.

# Version 1.0.4:
  * Fix issues with uninitialized memory if OLS coefficients cannot be computed.

# Version 1.0.3:
  * Fix an issue with exact fits and if the M-scale estimate is (close to) 0.

# Version 1.0.2:
  * Initial release.
