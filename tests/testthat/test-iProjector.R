test_that("G must be a numeric matrix",{
  data(nhanes)
  outcomes <- outcomes_tibble(nhanes)
  G <- rbind(
    norm = 1,
    BMXWT = outcomes$BMXWT
  )

  eta <- c(norm = 1, BMXWT = mean(nhanes$BMXWT))

  expect_no_error(iProjector(G, eta, v = empirical(outcomes)))
  ## G as a character matrix
  G["BMXWT", 2] <- "A"
  expect_error(
    iProjector(G, eta, v = empirical(outcomes)),
    regexp = "`G` must be a numeric matrix"
  )
  ## G is a list
  G <- list(
    norm = 1,
    BMXWT = outcomes$BMXWT
  )
  expect_error(
    iProjector(G, eta, v = empirical(outcomes)),
    regexp = "`G` must be a numeric matrix"
  )
})

test_that("all entries in `G` must be finite",{
  data(nhanes)
  outcomes <- outcomes_tibble(nhanes)
  G <- rbind(
    norm = 1,
    BMXWT = outcomes$BMXWT
  )

  eta <- c(norm = 1, BMXWT = mean(nhanes$BMXWT))

  expect_no_error(iProjector(G, eta, v = empirical(outcomes)))
  ## G contains NA
  G["BMXWT", 2] <- NA
  expect_error(
    iProjector(G, eta, v = empirical(outcomes)),
    regexp = "all entries in `G` must be finite"
  )
  ## G contains Inf
  G <- rbind(
    norm = Inf,
    BMXWT = outcomes$BMXWT
  )

  expect_error(
    iProjector(G, eta, v = empirical(outcomes)),
    regexp = "all entries in `G` must be finite"
  )
  ## G contains NaN
  G <- rbind(
    norm = NaN,
    BMXWT = outcomes$BMXWT
  )

  expect_error(
    iProjector(G, eta, v = empirical(outcomes)),
    regexp = "all entries in `G` must be finite"
  )
})

test_that("`eta` must be a numeric vector",{
  data(nhanes)
  outcomes <- outcomes_tibble(nhanes)
  G <- rbind(
    norm = 1,
    BMXWT = outcomes$BMXWT
  )

  eta <- c(norm = 1, BMXWT = mean(nhanes$BMXWT))

  expect_no_error(iProjector(G, eta, v = empirical(outcomes)))

  eta <- c(norm = "one", BMXWT = mean(nhanes$BMXWT))
  expect_error(
    iProjector(G, eta, v = empirical(outcomes)),
    regexp = "`eta` must be a numeric vector"
  )
  eta <- list(
    norm = 1,
    BMXWT = mean(nhanes$BMXWT)
  )
  expect_error(
    iProjector(G, eta, v = empirical(outcomes)),
    regexp = "`eta` must be a numeric vector"
  )
})

test_that("all entries in `eta` must be finite",{
  data(nhanes)
  outcomes <- outcomes_tibble(nhanes)
  G <- rbind(
    norm = 1,
    BMXWT = outcomes$BMXWT
  )

  eta <- c(norm = 1, BMXWT = mean(nhanes$BMXWT))

  expect_no_error(iProjector(G, eta, v = empirical(outcomes)))
  ## eta contains NA
  eta <- c(norm = NA, BMXWT = mean(nhanes$BMXWT))
  expect_error(
    iProjector(G, eta, v = empirical(outcomes)),
    regexp = "all entries in `eta` must be finite"
  )
  ## eta contains Inf
  eta <- c(norm = Inf, BMXWT = mean(nhanes$BMXWT))
  expect_error(
    iProjector(G, eta, v = empirical(outcomes)),
    regexp = "all entries in `eta` must be finite"
  )
  ## eta contain NaN
  eta <- c(norm = NaN, BMXWT = mean(nhanes$BMXWT))
  expect_error(
    iProjector(G, eta, v = empirical(outcomes)),
    regexp = "all entries in `eta` must be finite"
  )
})

test_that("`v` must be a numeric vector",{
  data(nhanes)
  outcomes <- outcomes_tibble(nhanes)
  G <- rbind(
    norm = 1,
    BMXWT = outcomes$BMXWT
  )
  eta <- c(norm = 1, BMXWT = mean(nhanes$BMXWT))
  expect_no_error(iProjector(G, eta, v = empirical(outcomes)))

  ## v is a character
  v <- c("A", 1)
  expect_error(
    iProjector(G, eta, v = v),
    regexp = "`v` must be a numeric vector"
  )
  ## v is a list
  v <- list(
    norm = 1,
    BMXWT = outcomes$BMXWT
  )
  expect_error(
    iProjector(G, eta, v = v),
    regexp = "`v` must be a numeric vector"
  )

})

test_that("all entries in `v` must be finite",{
  data(nhanes)
  outcomes <- outcomes_tibble(nhanes)
  G <- rbind(
    norm = 1,
    BMXWT = outcomes$BMXWT
  )
  eta <- c(norm = 1, BMXWT = mean(nhanes$BMXWT))
  v <- empirical(outcomes)
  expect_no_error(iProjector(G, eta, v = v))

  ## v contain NA
  v[1] <- NA
  expect_error(
    iProjector(G, eta, v = v),
    regexp = "all entries in `v` must be finite"
  )
  ## v contain Inf
  v[1] <- Inf
  expect_error(
    iProjector(G, eta, v = v),
    regexp = "all entries in `v` must be finite"
  )
  ## v contain NaN
  v[1] <- NaN
  expect_error(
    iProjector(G, eta, v = v),
    regexp = "all entries in `v` must be finite"
  )
})

test_that("`v` must be non-negative",{
  data(nhanes)
  outcomes <- outcomes_tibble(nhanes)
  G <- rbind(
    norm = 1,
    BMXWT = outcomes$BMXWT
  )
  eta <- c(norm = 1, BMXWT = mean(nhanes$BMXWT))
  v <- empirical(outcomes)
  expect_no_error(iProjector(G, eta, v = v))

  ## v has a negative entry
  v[1] <- -1
  expect_error(
    iProjector(G, eta, v = v),
    regexp = "`v` must be non-negative"
  )
})

test_that("`v` should sum to 1. Normalizing to 1",{
  data(nhanes)
  outcomes <- outcomes_tibble(nhanes)
  G <- rbind(
    norm = 1,
    BMXWT = outcomes$BMXWT
  )
  eta <- c(norm = 1, BMXWT = mean(nhanes$BMXWT))
  v <- empirical(outcomes)
  expect_no_error(iProjector(G, eta, v = v))

  ## v does not sum up to one
  v[1] <- 0
  expect_warning(
    iProjector(G, eta, v = v),
    regexp = "`v` should sum to 1. Normalizing to 1"
  )
})

test_that("`eta` must have one value per row of `G`", {
  data(nhanes)
  outcomes <- outcomes_tibble(nhanes)
  G <- rbind(
    norm = 1,
    BMXWT = outcomes$BMXWT
  )
  eta <- c(norm = 1, BMXWT = mean(nhanes$BMXWT))
  v <- empirical(outcomes)
  expect_no_error(iProjector(G, eta, v = v))

  ## eta has more entries than G has rows
  eta <- c(norm = 1, BMXWT = mean(nhanes$BMXWT), BMXHT = mean(nhanes$BMXHT))
  expect_error(
    iProjector(G, eta, v = v),
    regexp = "`eta` must have one value per row of `G`"
  )
  ## eta has less entries than G has rows
  eta <- c(norm = 1)
  expect_error(
    iProjector(G, eta, v = v),
    regexp = "`eta` must have one value per row of `G`"
  )
  ## G has more rows than eta has entries
  eta <- c(norm = 1, BMXWT = mean(nhanes$BMXWT))
  G <- rbind(
    norm = 1,
    BMXWT = outcomes$BMXWT,
    BMXHT = outcomes$BMXHT
  )
  expect_error(
    iProjector(G, eta, v = v),
    regexp = "`eta` must have one value per row of `G`"
  )
  ## G has less rows than eta has entries
  G <- G[1, ]
  expect_error(
    iProjector(G, eta, v = v),
    regexp = "`eta` must have one value per row of `G`"
  )

})
test_that("`v` must have one value per column of `G`",{
  data(nhanes)
  outcomes <- outcomes_tibble(nhanes)
  G <- rbind(
    norm = 1,
    BMXWT = outcomes$BMXWT
  )
  eta <- c(norm = 1, BMXWT = mean(nhanes$BMXWT))
  v <- empirical(outcomes)
  expect_no_error(iProjector(G, eta, v = v))

  ## v has less entries than G columns
  expect_error(
    iProjector(G, eta, v = v[-1] / sum(v[-1])),
    regexp = "`v` must have one value per column of `G`"
  )

  ## v has more entries than G columns
  expect_error(
    iProjector(G, eta, v = c(v, 1) / (sum(v) + 1)),
    regexp = "`v` must have one value per column of `G`"
  )
})

test_that("normalization condition is missing",{
  data(nhanes)
  outcomes <- outcomes_tibble(nhanes)
  G <- rbind(
    norm = 1,
    BMXWT = outcomes$BMXWT
  )
  eta <- c(norm = 1, BMXWT = mean(nhanes$BMXWT))
  v <- empirical(outcomes)
  expect_no_error(iProjector(G, eta, v = v))
  ## G misses normalization
  expect_error(
    iProjector(G[-1, , drop = FALSE], eta = eta[-2], v = v),
    regexp = "normalization condition is missing"
  )
  ## eta misses normalization
  expect_error(
    iProjector(G[-2, , drop = FALSE], eta = eta[-1], v = v),
    regexp = "normalization condition is missing"
  )
})

test_that("`G` has not full row rank",{
  data(nhanes)
  outcomes <- outcomes_tibble(nhanes)
  G <- rbind(
    norm = 1,
    BMXWT = outcomes$BMXWT
  )
  eta <- c(norm = 1, BMXWT = mean(nhanes$BMXWT))
  v <- empirical(outcomes)
  expect_no_error(iProjector(G, eta, v = v))
  ## one row is duplicated
  G <- rbind(
    norm = 1,
    BMXWT = outcomes$BMXWT,
    BMXWT2 = outcomes$BMXWT
  )
  eta <- c(norm = 1, BMXWT = mean(nhanes$BMXWT), BMXWT = mean(nhanes$BMXWT))
  expect_error(
    iProjector(G, eta = eta, v = v),
    regexp = "`G` has not full row rank"
  )
  ## normalization plus two marginals
  G <- rbind(
    norm = 1,
    male = outcomes$RIAGENDR == "male",
    female = outcomes$RIAGENDR == "female"
  )
  eta <- c(norm = 1, male = 0.5, female = 0.5)
  expect_error(
    iProjector(G, eta = eta, v = v),
    regexp = "`G` has not full row rank"
  )

})

test_that("`eta` is not a feasible",{
  data(nhanes)
  outcomes <- outcomes_tibble(nhanes)
  G <- rbind(
    norm = 1,
    BMXWT = outcomes$BMXWT
  )
  eta <- c(norm = 1, BMXWT = mean(nhanes$BMXWT))
  v <- empirical(outcomes)
  expect_no_error(iProjector(G, eta, v = v))
  ## mean weight of zero is not feasible
  eta = c(norm = 1, BMXWT = 0)
  expect_warning(
    iProjector(G, eta = eta, v = v),
    regexp = "`eta` is not feasible"
  )
  iProj <- suppressWarnings(iProjector(G, eta = eta, v = v))
  ## if eta is not feasible it returns v as I-projection
  expect_equal(iProj$p, v)
  ## converged is zero
  expect_equal(iProj$converged, 0)
  ## error is `eta` is not feasible
  expect_equal(iProj$error, "`eta` is not feasible")
  ## mean weight of max(BMXWT) + 0.1
  eta = c(norm = 1, BMXWT = max(outcomes$BMXWT) + 0.1)
  expect_warning(
    iProjector(G, eta = eta, v = v),
    regexp = "`eta` is not feasible"
  )
  ## mean weight of min(BMXWT) - 0.1
  eta = c(norm = 1, BMXWT = min(outcomes$BMXWT) - 0.1)
  expect_warning(
    iProjector(G, eta = eta, v = v),
    regexp = "`eta` is not feasible"
  )
  ## mean weight of max(BMXWT) is feasible --> no warning
  eta = c(norm = 1, BMXWT = max(outcomes$BMXWT))
  expect_no_warning(
    iProjector(G, eta = eta, v = v)
  )
  ## mean weight of min(BMXWT) is feasible --> no warning
  eta = c(norm = 1, BMXWT = min(outcomes$BMXWT))
  expect_no_warning(
    iProjector(G, eta = eta, v = v)
  )

})
