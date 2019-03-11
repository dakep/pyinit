#' Robust M-estimate of Scale
#'
#' Compute the M-estimate of scale using the MAD as initial estimate.
#'
#' This solves the M-estimation equation given by
#' \deqn{\sum_{i=1}^n \rho( x_i / s_n; cc ) = n delta}
#'
#' All \code{NA} values in \code{x} are removed before calculating the scale.
#'
#' @param x numeric vector.
#' @param delta desired value for the right-hand side of the M-estimation equation.
#' @param rho rho function to use in the M-estimation equation. Valid options
#' are \code{bisquare}, \code{huber} and \code{gauss}.
#' @param cc non-negative constant for the chosen rho function. If missing, it will be
#'          chosen such that the expected value of the rho function under the normal model
#'          is equal to \code{delta}.
#' @param eps threshold for convergence. Defaults to \code{1e-8}.
#' @param maxit maximum number of iterations. Defaults to \code{200}.
#'
#' @return Numeric vector of length one containing the solution \code{s_n} to
#' the equation above.
#'
#' @useDynLib pyinit C_mscale
mscale <- function(x, delta = 0.5, rho = c("bisquare", "huber", "gauss"), cc,
                   eps = 1e-8, maxit = 200) {

    if (!is.numeric(x) || !is.null(dim(x)) || length(x) == 0) {
        stop("`x` must be a numeric vector.")
    }
    if (any(is.na(x))) {
        x <- x[!is.na(x)]
    }

    rho <- match.arg(rho)

    if (missing(cc)) {
        cc <- NULL
    }

    ctrl <- pyinit_control(
        maxit = 1L,
        eps = 1e-6,
        psc_keep = 0.5,
        resid_keep_method = "threshold",
        resid_keep_prop = 2,
        resid_keep_thresh = 2,
        delta = delta,
        cc = cc,
        mscale_maxit = maxit,
        mscale_tol = eps,
        mscale_rho_fun = rho
    )

    x <- as.numeric(x)

    scale <- .Call(
        C_mscale,
        x,
        ctrl$delta,
        ctrl$mscale_cc,
        ctrl$mscale_maxit,
        ctrl$mscale_tol,
        ctrl$mscale_rho_fun
    )

    if (!is.finite(scale)) {
        scale <- NA_real_
    }

    return(scale)
}
