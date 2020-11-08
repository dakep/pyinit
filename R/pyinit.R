#' PY (Pena-Yohai) initial estimates for S-estimates of regression
#'
#' Computes the PY initial estimates for S-estimates of regression.
#'
#' @param x a matrix with the data, each observation in a row.
#' @param y the response vector.
#' @param intercept logical, should an intercept be included in the model? Defaults to
#'      \code{TRUE}.
#' @param delta,cc parameters for the M-scale estimator equation. If \code{cc} is
#'      missing it will be set to yield consistency under the Normal model for the
#'      given \code{delta} (right-hand side of the M-scale equation).
#' @param psc_keep proportion of observations to keep based on PSCs.
#' @param resid_keep_method how to clean the data based on large residuals.
#'      If \code{"threshold"}, all observations with scaled residuals larger
#'      than \code{resid_keep_thresh} will be removed (\code{resid_keep_thresh}
#'      corresponds to the constant \eqn{C_1} from equation (21) in Pena & Yohai
#'      (1999). If \code{"proportion"}, observations with the largest
#'      \code{resid_keep_prop} residuals will be removed.
#' @param maxit the maximum number of iterations to perform.
#' @param eps the relative tolerance for convergence. Defaults to \code{1e-8}.
#' @param resid_keep_prop,resid_keep_thresh see parameter
#'      \code{resid_keep_method} for details.
#' @param mscale_maxit maximum number of iterations allowed for the M-scale
#'      algorithm. Defaults to \code{200}.
#' @param mscale_tol convergence threshold for the m-scale
#' @param mscale_rho_fun A string containing the name of the rho
#' function to use for the M-scale. Valid options
#' are \code{bisquare}, \code{huber} and \code{gauss}.
#'
#' @return
#' \item{coefficients}{numeric matrix with coefficient vectors in columns. These
#' are regression estimators based on "cleaned" subsets of the data. The
#' M-scales of the corresponding residuals are returned in the entry
#' \code{objective}. The regression coefficients with smallest estimated
#' residual scale is in the first column, but the others need not be ordered.}
#' \item{objective}{vector of values of the M-scale estimate of the residuals
#' associated with each vector of regression coefficients in the columns of
#' \code{coefficients}.}
#'
#' @references Pena, D., & Yohai, V.. (1999). A Fast Procedure for Outlier Diagnostics in Large
#' Regression Problems. \emph{Journal of the American Statistical Association}, 94(446),
#' 434-445. <doi:10.2307/2670164>
#'
#' @examples
#' # generate a simple synthetic data set for a linear regression model
#' # with true regression coefficients all equal to one "(1, 1, 1, 1, 1)"
#' set.seed(123)
#' x <- matrix(rnorm(100*4), 100, 4)
#' y <- rnorm(100) + rowSums(x) + 1
#' # add masked outliers
#' a <- svd(var(x))$v[,4]
#' x <- rbind(x, t(outer(a, rnorm(20, mean=4, sd=1))))
#' y <- c(y, rnorm(20, mean=-2, sd=.2))
#'
#' # these outliers are difficult to find
#' plot(lm(y~x), ask=FALSE)
#'
#' # use pyinit to obtain estimated regression coefficients
#' tmp <- pyinit(x=x, y=y, resid_keep_method='proportion', psc_keep = .5, resid_keep_prop=.5)
#' # the vector of regression coefficients with smallest residuals scale
#' # is returned in the first column of the "coefficients" element
#' tmp$coefficients[,1]
#' # compare that with the LS estimator on the clean data
#' coef(lm(y~x, subset=1:100))
#' # compare it with the LS estimator on the full data
#' coef(lm(y~x))
#'
#'
#' @useDynLib pyinit C_initpy
#' @importFrom stats var
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
    mscale_rho_fun = c("bisquare", "huber", "gauss")
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
            coefficients = matrix(est, nrow = 1L, ncol = 1L),
            objective = mscale(
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
        coefficients = initCoefs,
        objective = ies[[3L]][seq_len(ies[[1L]])]
    ))
}
