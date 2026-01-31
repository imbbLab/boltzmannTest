#' @export
boltzmannTest <- function(x, ...) UseMethod("boltzmannTest")

## what do we need
## 1. G the coefficient matrix
## 2. H = eta == vector of hypothesis conditions
## 3. A = mu == vector of observed means
## 4. A' = subset of observed means, confounders?

## boltzmannTest <- function(G, eta, mu, subs)
## 1. project empirical to hypothesis linear family H: G p = eta -> h
## 2. if subs are present project h to to A_k: G_k p = mu_k -> a_k
## 3. if subs are not present a_k = h (in a sense only normalization)
## 4. project a_k to A: G * p = mu -> f

#' @export
boltzmannTest.entities_tibble <-function(entities, G, eta, nested = NULL, strict = FALSE, maxit = 10000L, tolerance = .Machine$double.eps){
  if (NCOL(G) != length(eta)){
    stop("the number of constraints ")
  }
  if (!is.null(nested)){
    nested <- as.vector(nested)
    if (any(nested) < 1 || any(nested) > NROW(G)){
      stop("the indices of the nested conditions are out of range")
    }
  }
  ## determine the generalized moments for the empirical distribution
  mu <- G %*% empirical(entities)

  ## check for equality per entry
  delta <- mu - eta
  ## need approximate equality, so we are testing the absolute value
  ## of each entry of the difference vector delta
  isDifferent <- sapply(
    delta,
    function(x){
      ! abs(x) < tolerance
    }
  )
  ## if they are different they constitute the to be tested generalized moments, i.e. their number
  ## corresponds to the degrees of freedom of the chi^2
  dof <- sum(isDifferent)
  ## if they are equal they probably constitute structural conditions (like normalization)
  ## or some fuckup on the side of the user
  ##
  ## for nested testing we need further information about the "different" generalized moments
  ## that first need to be sampled from the hypothesis distribution...


  ## 1. project empirical distribution to hypothesis linear family
  h <- iProjector(G = G, eta = eta, v = empirical(entities), maxit = maxit, convTolerance = tolerance)
  ## check for convergence
  if(h$converged != 0){
    ## "update" hypothesis distribution by the empirical moments
    if (!is.null(nested)){

    }
    ## here we cover the case, when some p's got zero
    a <- if(any(h$p < tolerance)){
      stop("hypothesis conditions led to zero probabilities")
      ## we may remove the offending entry and try again?
    }
    idiv <- iDivergence(empirical(entities), h$p)

    statistic <- 2 * sampleSize(entities) * idiv
    p.value = pchisq(statistic, df = dof, lower.tail = FALSE)




  }




}

res = list(
  statistic = 10,
  parameter = c(idiv = 0.1, dof = 2L),
  p.value = 0.03,
  estimate = c(0.1, 0.2),
  null.value = c(0.0, 0.0),
  method = "Boltzmann Test",
  data.name = "data"
)
class(res) = "htest"
