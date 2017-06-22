## Creates the internal control list for PY initial estimators
##
## Takes care of the correct storage mode for the arguments passed to the C/C++ code
##
## @param lambda numeric >= 0
## @param alpha numeric 0 <= alpha <= 1
## @param numIt integer > 0
## @param eps numeric > 0
## @param resid.clean.method character
## @param resid.threshold If \code{resid.clean.method = "threshold"} numeric > 0, otherwise not
##          referenced
## @param resid.proportion If \code{resid.clean.method = "proportion"} numeric > 0, otherwise not
##          referenced
## @param psc.proportion numeric > 0
## @param mscale.delta numeric > 0
## @param mscale.cc numeric > 0
initest.control <- function(
    numIt,
    eps = 1e-6,
    resid.clean.method = c("proportion", "threshold"),
    resid.threshold,
    resid.proportion,
    psc.proportion,
    mscale.delta,
    mscale.cc,
    mscale.maxit,
    mscale.tol,
    mscale.rho.fun
) {
    ret <- as.list(environment())

    simpleCheck <- function(x, min = 0, max = Inf, eqMin = FALSE, eqMax = FALSE) {
        if (length(x) != 1L || !is.numeric(x) || is.na(x) || x < min || x > max ||
            (!eqMin && (x == min)) || (!eqMax && (x == max))) {
            stop(sprintf("`%s` must be single number between %f and %f",
                         deparse(substitute(x)), min, max))
        }
    }

    ret$resid.clean.method <- match.arg(resid.clean.method)
    ret$mscale.rho.fun <- .rho2IntRho(ret$mscale.rho.fun)

    if (ret$resid.clean.method == "proportion") {
        simpleCheck(resid.proportion)

        ret$resid.threshold <- -1

        if (resid.proportion > 1) {
            stop("`resid.proportion` must be less than 1")
        }
    } else {
        simpleCheck(resid.threshold)
        ret$resid.proportion <- -1
    }

    simpleCheck(numIt)
    simpleCheck(eps)
    simpleCheck(psc.proportion)
    simpleCheck(mscale.delta)

    if (mscale.cc <= 0) {
        ret$mscale.cc <- consistency.rho(ret$mscale.delta, ret$mscale.rho.fun)
    }

    simpleCheck(ret$mscale.cc)
    simpleCheck(mscale.maxit)
    simpleCheck(mscale.tol)
    simpleCheck(ret$mscale.rho.fun, eqMin = TRUE)

    if (ret$psc.proportion > 1) {
        stop("`psc.proportion` must be less than 1")
    }

    ret$numIt <- as.integer(ret$numIt)
    ret$mscale.maxit <- as.integer(ret$mscale.maxit)

    return(ret)
}

.rho2IntRho <- function(rho.fun) {
    rho.fun <- as.integer(pmatch(rho.fun, c("bisquare", "huber"))) - 1L

    rho.fun <- rho.fun[which(!is.na(rho.fun))[1L]]

    if (is.na(rho.fun)) {
        rho.fun <- 0L
        warning("Unknown rho function selected. Using Tukey's bisquare.")
    }

    return(rho.fun)
}
