library(pyinit)
library(testthat)

test_that("Exact-fit problems", {
    set.seed(123)
    x <- rnorm(50)
    y <- 2 + 3 * x
    X <- cbind(1, x)

    n <- nrow(X)
    p <- ncol(X)

    res <- pyinit(x = X, y = y, intercept = FALSE, delta = 0.5 * (1 - p/n), cc = 1.54764,
                  psc_keep = 0.5 * (1 - p/n), resid_keep_method = 'threshold',
                  resid_keep_thresh = 2, resid_keep_prop = 0.2, maxit = 20L,  eps = 1e-5,
                  mscale_maxit = 50, mscale_tol = 1e-6, mscale_rho_fun = 'bisquare')

    expect_equal(c(0, 0), res$objective, tolerance = 1e-5,
                 info = 'Objective values are different')
    expect_equal(cbind(c(2, 3), c(2, 3)), res$coefficients, tolerance = 1e-5,
                 info = 'Coefficients are different')
})


test_that("Rank-deficient problems", {
    X <- cbind(1, mtcars$wt, mtcars$gear)
    y <- as.vector(mtcars$mpg)
    n <- nrow(X)
    p <- ncol(X)

    res <- pyinit(x = X, y = y, intercept = FALSE, delta = 0.5 * (1 - p/n), cc = 1.54764,
                  psc_keep = 0.5 * (1 - p/n), resid_keep_method = 'threshold',
                  resid_keep_thresh = 2, resid_keep_prop = 0.2, maxit = 20L,  eps = 1e-5,
                  mscale_maxit = 50, mscale_tol = 1e-6, mscale_rho_fun = 'bisquare')

    expected_obj <- c(3.16984, 3.16984, 3.22770, 3.24876, 3.40174, 3.43364, 3.47512, 3.66026, 5.07470, 5.07470)
    obj_order <- order(res$objective)
    expect_equal(expected_obj, res$objective[obj_order], tolerance = 1e-5,
                 info = 'Objective values are different')

    expect_known_hash(round(res$coefficients[, obj_order], 4), '30f3b173bb32999ace3f3072ed')
})
