#' PY (Pena-Yohai) initial estimates for S-estimates of regression
#'
#' Computes the PY initial estimates for S-estimates of regression.
#'
#' @param x the data matrix X.
#' @param y the response vector.
#' @param intercept should an intercept be included in the models. Defaults to \code{TRUE}.
#' @param delta,cc parameters for the M-equation of the scale. If \code{cc} is
#'      missing it will be set to yield consistency under the Normal model.
#' @param psc_keep the proportion of observations to keep based on PSCs.
#' @param resid_keep_method how to clean the data based on large residuals.
#'      If \code{"threshold"}, all observations with scaled residuals larger than
#'      \code{resid_keep_thresh} will be removed (\code{resid_keep_thresh} corresponds to the constant
#'      \eqn{C_1} from equation (21) in Pena & Yohai (1999).
#'      If \code{"proportion"}, observations with the largest \code{resid_keep_prop} residuals
#'      will be removed.
#' @param maxit the maximum number of iterations to perform.
#' @param eps the relative tolerance for convergence.
#' @param resid_keep_prop,resid_keep_thresh see parameter \code{resid_keep_method} for details.
#' @param mscale_maxit Maximum number of iterations allowed for the M-scale algorithm.
#' @param mscale_tol Convergence threshold for the m-scale
#' @param mscale_rho_fun The rho function to use for the m-scale.
#'
#' @return \item{initCoef}{A numeric matrix with one initial coefficient per column}
#'         \item{objF}{A vector of values of the objective function for the respective coefficient}
#'
#' @references Pena, D., & Yohai, V.. (1999). A Fast Procedure for Outlier Diagnostics in Large
#' Regression Problems. \emph{Journal of the American Statistical Association}, 94(446),
#' 434-445. \url{http://doi.org/10.2307/2670164}
#'
#' @useDynLib pyinit C_initpy
#' @export
pyinit <- function(
    x, y,
    intercept = TRUE,
    delta = 0.5,
    cc,
    maxit = 10,
    psc_keep,
    resid_keep_method = c("threshold", "proportion"),
    resid_keep_prop,
    resid_keep_thresh,
    eps = 1e-8,
    mscale_maxit = 200,
    mscale_tol = eps,
    mscale_rho_fun = c("bisquare", "huber")
) {
    y <- drop(y)

    dx <- dim(x)
    dY <- dim(y)
    yl <- length(y)

    if (is.null(yl) || (!is.null(dY) && length(dY) != 1L) || !is.numeric(y)) {
        stop("`yl` must be a numeric vector")
    }

    # Force y to be stored as double
    storage.mode(y) <- "double"

    if (is.null(dx) || length(dx) != 2L || !is.numeric(x) || dx[1L] != yl) {
        stop("`x` must be a numeric matrix with the same number of observations as `y`")
    }

    if (anyNA(x) || anyNA(y)) {
        stop("Missing values are not supported")
    }

    if (!identical(intercept, FALSE)) {
        ## Add leading column of 1's
        x <- cbind(1, x)
        dx[2L] <- dx[2L] + 1L
    }

    # Force x to be stored as double
    storage.mode(x) <- "double"
    resid_keep_method <- match.arg(resid_keep_method)

    if (resid_keep_method == "threshold") {
        resid_keep_prop <- 0
    } else {
        resid_keep_thresh <- 0
    }

    if (missing(cc)) {
        cc <- NULL
    }

    mscale_rho_fun <- match.arg(mscale_rho_fun)

    ctrl <- pyinit_control(
        maxit = maxit,
        eps = eps,
        delta = delta,
        cc = cc,
        psc_keep = psc_keep,
        resid_keep_method = resid_keep_method,
        resid_keep_prop = resid_keep_prop,
        resid_keep_thresh = resid_keep_thresh,
        mscale_maxit = mscale_maxit,
        mscale_tol = mscale_tol,
        mscale_rho_fun = mscale_rho_fun
    )

    # Check if x only has a single constant column and compute the only
    # candidate by hand
    if ((dx[2L] == 1L) && (var(x) < .Machine$double.eps)) {
        est <- mean(y) / mean(x)
        return(list(
            initCoef = matrix(est, nrow = 1L, ncol = 1L),
            objF = mscale(
                y - est,
                delta = ctrl$delta,
                rho = mscale_rho_fun,
                cc = ctrl$mscale_cc,
                eps = ctrl$mscale_tol,
                maxit = ctrl$mscale_maxit
            )
        ))
    }

    usableProp <- 1
    if ((ctrl$resid_keep_method == "proportion") && (ctrl$maxit > 1L)) {
        usableProp <- ctrl$psc_keep * ctrl$resid_keep_prop
    }

    if (dx[2L] >= dx[1L]) {
        stop(
            "`pyinit` can not be used for data with more variables than ",
            "observations"
        )
    } else if (dx[2L] >= ceiling(usableProp * dx[1L])) {
        stop(
            "With the specified proportion of observations to remove, the ",
            "number of observations will be smaller than the number of ",
            "variables.\nIn this case `pyinit` can not be used."
        )
    }

    dx <- dim(x)

    ies <- .Call(
        C_initpy,
        t(x),
        y,
        ctrl$maxit,
        ctrl$eps,
        ctrl$resid_keep_thresh,
        ctrl$resid_keep_prop,
        ctrl$psc_keep,
        ctrl$delta,
        ctrl$mscale_cc,
        ctrl$mscale_maxit,
        ctrl$mscale_tol,
        ctrl$mscale_rho_fun
    )

    initCoefs <- matrix(
        ies[[2L]][seq_len(ies[[1L]] * dx[2L])],
        nrow = dx[2L], ncol = ies[[1L]]
    )

    return(list(
        initCoef = initCoefs,
        objF = ies[[3L]][seq_len(ies[[1L]])]
    ))
}
