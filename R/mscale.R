#' Robust M-estimate of Scale
#'
#' Compute the M-estimate of scale with MAD as initial estimate.
#'
#' This solves the M-estimation equation given by
#' \deqn{\sum_{i=1}^n \rho( x_i / s_n; cc ) = n b}
#'
#' \code{NA} values in \code{x} are removed before calculating the
#'
#' @param x A numeric vector.
#' @param delta The value of the M-estimation equation.
#' @param rho The rho function to use in the M-estimation equation.
#' @param cc Non-negative constant for the chosen rho function. If missing, it will be
#'          chosen such that the expected value of the rho function under the normal model
#'          is equal to \code{b}.
#' @param eps Threshold for convergence
#' @param max.it The maximum number of iterations
#'
#' @return Numeric vector of length one
#'
#' @useDynLib pyinit C_mscale
mscale <- function(x, delta = 0.5, rho = c("bisquare", "huber"), cc,
                   eps = 1e-8, max.it = 200) {

    if (!is.numeric(x) || !is.null(dim(x)) || length(x) == 0) {
        stop("`x` must be a numeric vector.")
    }
    if (any(is.na(x))) {
        x <- x[!is.na(x)]
    }

    rho <- match.arg(rho)

    if (missing(cc)) {
        cc <- 0
    }

    ctrl <- initest.control(
        numIt = 1L,
        eps = 1e-9,
        resid.clean.method = "proportion",
        resid.threshold = 2,
        resid.proportion = .5,
        psc.proportion = .5,
        mscale.delta = delta,
        mscale.cc = cc,
        mscale.maxit = max.it,
        mscale.tol = eps,
        mscale.rho.fun = rho
    )

    x <- as.numeric(x)

    scale <- .Call(
        C_mscale,
        x,
        length(x),
        ctrl$mscale.delta,
        ctrl$mscale.cc,
        ctrl$mscale.maxit,
        ctrl$mscale.tol,
        ctrl$mscale.rho.fun
    )

    if (!is.finite(scale)) {
        scale <- NA_real_
    }

    return(scale)
}
