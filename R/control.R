## Creates the internal control list for PY initial estimators
##
## Takes care of the correct storage mode for the arguments passed to the C/C++ code
##
pyinit_control <- function(
    maxit,
    eps,
    psc_keep,
    resid_keep_method = c("proportion", "threshold"),
    resid_keep_prop,
    resid_keep_thresh,
    delta,
    cc,
    mscale_maxit,
    mscale_tol,
    mscale_rho_fun = c("bisquare", "huber")
) {
    mscale_rho_fun <- .rho_to_int_rho(mscale_rho_fun)
    resid_keep_method <- match.arg(resid_keep_method)

    if (resid_keep_method == "threshold") {
        resid_keep_thresh <- .check_arg(resid_keep_thresh, "numeric", range = 0)
        resid_keep_prop <- -1
    } else {
        resid_keep_prop <- .check_arg(resid_keep_prop, "numeric", range = c(0, 1))
        resid_keep_thresh <- -1
    }

    delta <- .check_arg(delta, "numeric", range = c(0, 0.5), range_test_upper = "<=")
    if (missing(cc) || is.null(cc)) {
        cc <- .consistency_rho(delta, mscale_rho_fun)
    }

    return(list(
        maxit = .check_arg(maxit, "integer", range = 0),
        eps = .check_arg(eps, "numeric", range = 0),
        psc_keep = .check_arg(psc_keep, "numeric", range = c(0, 1)),
        resid_keep_method = resid_keep_method,
        resid_keep_prop = resid_keep_prop,
        resid_keep_thresh = resid_keep_thresh,
        delta = delta,
        mscale_cc = cc,
        mscale_maxit = .check_arg(mscale_maxit, "integer", range = 0),
        mscale_tol = .check_arg(mscale_tol, "numeric", range = 0),
        mscale_rho_fun = mscale_rho_fun
    ))
}



## Utility function to check arguments for the correct type and optionally
## length and value range.
##
##
.check_arg <- function (
    x,
    type = c("integer", "numeric", "logical", "character"),
    range,
    length = 1L,
    range_test_lower = ">",
    range_test_upper = "<"
) {
    type <- match.arg(type)
    ok_type <- switch(
        type,
        integer = is.integer(x) || (
            is.numeric(x) && isTRUE(abs(as.integer(x) - x) < .Machine$double.eps)
        ),
        numeric = is.numeric(x),
        logical = is.logical(x),
        character = is.character(x)
    )

    if (!isTRUE(ok_type)) {
        stop(sprintf("%s is expected to be `%s`", deparse(substitute(x)), type))
    }

    x <- switch(
        type,
        integer = as.integer(x),
        numeric = as.double(x),
        x
    )

    ok_length <- if (!is.null(length)) {
        length(x) == length
    } else {
        TRUE
    }

    if (!isTRUE(ok_length)) {
        stop(sprintf("%s is expected to be of length %d",
                     deparse(substitute(x)), length))
    }

    range_fun_lower <- match.fun(range_test_lower)
    range_fun_upper <- match.fun(range_test_upper)

    ok_range <- if (!missing(range)) {
        if (length(range) == 1L) {
            all(range_fun_lower(x, range[1L]))
        } else {
            all(range_fun_lower(x, range[1L]) & range_fun_upper(x, range[2L]))
        }
    } else {
        TRUE
    }

    if (!isTRUE(ok_range)) {
        message <- if (length(range) > 1L) {
            message <- paste(
                range_test_lower,
                range[1L],
                "and",
                range_test_upper,
                range[2L],
                sep = " "
            )
        } else {
            paste(range_test_lower, range[1L], sep = " ")
        }
        stop(sprintf("%s is expected to be %s",
                     deparse(substitute(x)), message))
    }

    return(x)
}


## Get the constant needed for consistency for the given delta
## and the given rho function
#' @importFrom robustbase .Mchi
#' @importFrom stats dnorm pnorm integrate uniroot
.consistency_rho <- function(delta, int_rho_fun) {
    if (is.character(int_rho_fun)) {
        int_rho_fun <- .rho_to_int_rho(int_rho_fun)
    }

    ##
    ## Pre-computed values for some delta values
    ##
    if (abs(delta - 0.5) < sqrt(.Machine$double.eps)) {
        return(switch(
            int_rho_fun + 1L,
            1.368482, # huber
            1.547645  # bisquare
        ))
    } else if (abs(delta - 0.25) < sqrt(.Machine$double.eps)) {
        return(switch(
            int_rho_fun + 1L,
            1.988013, # huber
            2.937015  # bisquare
        ))
    } else if (abs(delta - 0.1) < sqrt(.Machine$double.eps)) {
        return(switch(
            int_rho_fun + 1L,
            3.161931, # huber
            5.182361  # bisquare
        ))
    } else if (delta < 0.005) {
        return(50) # ~.1% bdp for bisquare, 9.6e-5% for huber
    }

    integral_interval <- c(1, 7)

    if (delta < 0.1) {
        integral_interval <- c(3, 30)
    }

    integrand_huber <- function(x, cc) {
        dnorm(x) * .Mchi(x, cc, 1L) / (0.5 * cc * cc)
    }

    expectation <- if (int_rho_fun == 1L) {
        # For bisquare we have the closed form solution to the expectation
        function(cc, delta) {
            pnorm.mcc <- 2 * pnorm(-cc)
            1/cc^6 * exp(-(cc^2/2)) * (
                -cc * (15 - 4 * cc^2 + cc^4) * sqrt(2 / pi) +
                    3 * (5 - 3 * cc^2 + cc^4) * exp(cc^2/2) * (1 - pnorm.mcc) +
                    cc^6 * exp(cc^2/2) * pnorm.mcc
            ) - delta
        }
    } else {
        function(cc, delta) {
            integrate(integrand_huber, lower = -Inf, upper = Inf, cc)$value - delta
        }
    }

    uniroot(expectation, interval = integral_interval, delta)$root
}


.rho_to_int_rho <- function(rho_fun = c("huber", "bisquare")) {
    rho_fun <- switch(
        match.arg(rho_fun),
        huber = 0L,
        bisquare = 1L
    )

    rho_fun <- rho_fun[which(!is.na(rho_fun))[1L]]

    if (is.na(rho_fun)) {
        rho_fun <- 1L
        warning("Unknown rho function selected. Using Tukey's bisquare.")
    }

    return(rho_fun)
}
