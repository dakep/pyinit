
## Get the constant needed for consistency for the given delta
## and the given rho function
#' @importFrom robustbase .Mchi
consistency.rho <- function(delta, int.rho.fun, interval = c(0.3, 10)) {
    if (is.character(int.rho.fun)) {
        int.rho.fun <- switch (int.rho.fun,
                               huber = 0L,
                               bisquare = 1L,
                               1L
        )
    }

    integrand <- function(x, cc) {
        dnorm(x) * .Mchi(x, cc, int.rho.fun)
    }

    expectation <- function(cc, delta) {
        integrate(integrand, lower = -Inf, upper = Inf, cc)$value - delta
    }

    uniroot(expectation, interval = interval, delta)$root
}
