#' @export
boltzmannTest <- function(x, ...) UseMethod("boltzmannTest")

## one and two sample tests for the mean
#' @importFrom tibble tibble
#' @export
boltzmannTest.default <-function(x, y = NULL, mu = 0, paired = FALSE){
  if (!missing(mu) && (length(mu) != 1 || us.na(mu))){
    stop("`mu` must be a single number")
  }
  if (!is.numeric(x)){
    stop("`x` must be a numeric vector")
  }
  xName <- deparse(substitute(x))
  if(!is.null(y)){
    if (!is.numeric(y)){
      stop("`y` must a be a numeric vector")
    }
    yName <- deparse(substitute(y))
    dataName <- paste0(xName, " and", yName)
    if (paired){
      if(length(x) != length(y)){
        stop(paste0("`", xName, "` and `", yName, "` must have the same length"))
      }
      method <- "Paired Boltzmann Test of a mean"

      ok <-complete.cases(x,y)
      x <- x[ok]
      y <- y[ok]
      data <- tibble(
        x = x,
        y = y
      )
      data <- data[complete.cases(data), ]
      entities <- entities_tibble(data)
      G <- rbind(
        1,
        entities$x - entities$y
      )
      eta <- c(1, mu)
      names(eta) = c("norm", paste0("<", xName, "-", yName, ">"))
      testedMoments <- 2
    } else{
      method <- "Two Sample Boltzmann Test of the difference of means"
      xok <- !is.na(x)
      yok <- !is.na(y)

      data <- tibble(
        z = c(x[xok], y[yok]),
        group = factor(c(rep("x", sum(xok)), rep("y", sum(yok))), levels = c("x", "y"))
      )
      entities <- entities_tibble(data)
      xNameProb <- sum((entities$group == "x") * empirical(entities))
      G <- rbind(
        1,
        with(
          entities,
          ifelse(
            group == "x",
            z / xNameProb,
            -z / (1 - xNameProb)
          )
        ),
        with(
          entities,
          group == "x"
        )
      )
      eta <- c(1, mu, xNameProb)
      names(eta) = c("norm", paste0("<x_", xName, ">-<y_", yName, ">"))
      testedMoments <- 2
    }
  } else{

    if (paired){
      stop("`y` is missing for paired test")
    }
    dataName <- xName
    method <- "One Sample Boltzmann Test of the mean"
    xok <- !is.na(x)
    x <- x[xok]
    data <- tibble(
      z = x
    )
    entities <- entities_tibble(data)

    G <- rbind(
      1,
      entities$z
    )

    eta <- c(1, mu)
    names(eta) <- c("norm", paste0("<", xName, ">"))
    testedMoments <- 2
  }
  res <- boltzmannTest(entities, G, eta, testedMoments)
  res$dataName = dataName
  res$method = method
  res
}



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
boltzmannTest.entities_tibble <-function(entities, G, eta, testedMoments, nested = NULL, strict = FALSE, maxit = 10000L, tolerance = .Machine$double.eps){
  if (NROW(G) != length(eta)){
    stop("the number of generalized moments given by `eta` does not match the number of rows in `G`")
  }
  if (!is.null(nested)){
    nested <- as.vector(nested)
    if (any(nested) < 1 || any(nested) > NROW(G)){
      stop("the indices of the nested moments are out of range")
    }
    if (any(nested %in% testedMoments)){
      stop("the nested moments must not include the tested moments")
    }
  }
  testedMoments <- as.integer(testedMoments)
  if (any(testedMoments < 1) || any(testedMoments > NROW(G))){
    stop("the indices of the tested moments are out of range")
  }
  if (length(testedMoments) < 1){
    stop("there must be at least one moment to be tested")
  }
  degreesOfFreedom <- length(testedMoments)
  dataName <- deparse(substitute(entities))
  sampleSize <- sampleSize(entities)
  ## determine the generalized moments for the empirical distribution
  mu <- G %*% empirical(entities)

  etaNested <- NULL
  if (!is.null(nested)){
    etaNested <- eta
    eta[nestedMoments] = mu[nestedMoments]
  }

  ## 1. project empirical distribution to hypothesis linear family
  h <- iProjector(G = G, eta = eta, v = empirical(entities), maxit = maxit, convTolerance = tolerance)
  ## check for convergence
  if(h$converged != 0){

    ## here we cover the case, when some p's got zero
    a <- if(any(h$p < tolerance)){
      stop("hypothesis conditions led to zero probabilities")
      ## we may remove the offending entry and try again?
    }
    ## "update" hypothesis distribution by the empirical moments

    baseDistribution <- if (!is.null(nested)){
      iProjector(G = G, eta = etaNested, v = h$p, maxit = maxit, convTolerance = tolerance)
    } else{
      h
    }
    iDivergence <- iDivergence(empirical(entities), baseDistribution$p)

    statistic <- 2 * sampleSize * iDivergence
    pValue = pchisq(statistic, df = degreesOfFreedom, lower.tail = FALSE)
    boltzmannTestResult(
      statistic = statistic,
      iDivergence = iDivergence,
      degreesOfFreedom = degreesOfFreedom,
      sampleSize = sampleSize,
      pValue = pValue,
      alternativeDistribution = empirical(entities),
      alternativeMoments = mu,
      hypothesisDistribution = h$p,
      hypothesisMoments = eta,
      testedMoments = testedMoments,
      dataName = dataName,
      coefficientMatrix = G,
      nestedAlternativeDistribution = if (is.null(nested)) NULL else baseDistribution$p,
      nestedAlternativeMoments = etaNested
    )
  } else{
    boltzmannTestResult(
      statistic = NA,
      iDivergence = NA,
      degreesOfFreedom = degreesOfFreedom,
      sampleSize = sampleSize,
      pValue = NA,
      alternativeDistribution = empirical(entities),
      alternativeMoments = mu,
      hypothesisDistribution = rep(NA, NROW(entities)),
      hypothesisMoments = eta,
      testedMoments = testedMoments,
      dataName = dataName,
      coefficientMatrix = G,
      nestedAlternativeDistribution = if (is.null(nested)) NULL else rep(NA, NROW(entities)),
      nestedAlternativeMoments = if (is.null(nested)) NULL else etaNested
    )
  }




}

