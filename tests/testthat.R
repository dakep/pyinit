if (require(testthat)) {
  library(pyinit)
  test_check("pyinit")
} else {
  warning("'pyinit' requires 'testthat' for tests.")
}
