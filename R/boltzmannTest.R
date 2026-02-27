#' @export
boltzmannTest <- function(x, ...) UseMethod("boltzmannTest")

#' Boltzmann Test
#'
#' @importFrom tibble tibble
#' @export
#'

boltzmannTest.outcomes_tibble <-function(outcomes, G, eta, testedMoments, nested = NULL, strict = FALSE, maxit = 10000L, tolerance = .Machine$double.eps){
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
  dataName <- deparse(substitute(outcomes))
  sampleSize <- sampleSize(outcomes)
  ## determine the generalized moments for the empirical distribution
  mu <- G %*% empirical(outcomes)

  etaNested <- NULL
  if (!is.null(nested)){
    etaNested <- eta
    eta[nestedMoments] = mu[nestedMoments]
  }

  ## 1. project empirical distribution to hypothesis linear family
  h <- iProjector(G = G, eta = eta, v = empirical(outcomes), maxit = maxit, convTolerance = tolerance)
  ## check for convergence
  if(h$converged != 0){

    ## here we cover the case, when some p's got zero
    a <- if(any(h$p < tolerance)){
      warning("hypothesis conditions led to zero probabilities")
      return(boltzmannTestResult(
        statistic = NA,
        iDivergence = NA,
        degreesOfFreedom = degreesOfFreedom,
        sampleSize = sampleSize,
        pValue = NA,
        alternativeDistribution = empirical(outcomes),
        alternativeMoments = mu,
        hypothesisDistribution = rep(NA, NROW(outcomes)),
        hypothesisMoments = eta,
        testedMoments = testedMoments,
        dataName = dataName,
        coefficientMatrix = G,
        nestedAlternativeDistribution = if (is.null(nested)) NULL else rep(NA, NROW(outcomes)),
        nestedAlternativeMoments = if (is.null(nested)) NULL else etaNested
      ))
      ## we may remove the offending entry and try again?
    }
    ## "update" hypothesis distribution by the empirical moments

    baseDistribution <- if (!is.null(nested)){
      iProjector(G = G, eta = etaNested, v = h$p, maxit = maxit, convTolerance = tolerance)
    } else{
      h
    }
    iDivergence <- iDivergence(empirical(outcomes), baseDistribution$p)

    statistic <- 2 * sampleSize * iDivergence
    pValue = pchisq(statistic, df = degreesOfFreedom, lower.tail = FALSE)
    boltzmannTestResult(
      statistic = statistic,
      iDivergence = iDivergence,
      degreesOfFreedom = degreesOfFreedom,
      sampleSize = sampleSize,
      pValue = pValue,
      alternativeDistribution = empirical(outcomes),
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
      alternativeDistribution = empirical(outcomes),
      alternativeMoments = mu,
      hypothesisDistribution = rep(NA, NROW(outcomes)),
      hypothesisMoments = eta,
      testedMoments = testedMoments,
      dataName = dataName,
      coefficientMatrix = G,
      nestedAlternativeDistribution = if (is.null(nested)) NULL else rep(NA, NROW(outcomes)),
      nestedAlternativeMoments = if (is.null(nested)) NULL else etaNested
    )
  }




}

#' @importFrom tibble tibble
#' @export
boltzmannTest.numeric <-function(x, y = NULL, mu = 0, paired = FALSE){
  if (!missing(mu) && (length(mu) != 1 || is.na(mu))){
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
      outcomes <- outcomes_tibble(data)
      G <- rbind(
        1,
        outcomes$x - outcomes$y
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
      outcomes <- outcomes_tibble(data)
      xNameProb <- sum((outcomes$group == "x") * empirical(outcomes))
      G <- rbind(
        1,
        with(
          outcomes,
          ifelse(
            group == "x",
            z / xNameProb,
            -z / (1 - xNameProb)
          )
        ),
        with(
          outcomes,
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
    outcomes <- outcomes_tibble(data)

    G <- rbind(
      1,
      outcomes$z
    )

    eta <- c(1, mu)
    names(eta) <- c("norm", paste0("<", xName, ">"))
    testedMoments <- 2
  }
  res <- boltzmannTest(outcomes, G, eta, testedMoments)
  res$dataName = dataName
  res$method = method
  res
}

#' @export
boltzmannTest.formula <- function(formula, data, eta = 0){
  if (missing(formula)){
    stop("`formula` not provided")
  }
  ## coerce data to a data frame
  data <- as.data.frame(data)

  ## make sure that all variables used in the formula are also present
  ## as columns in the data
  vars <- all.vars(formula)
  if (!all(vars %in% colnames(data))) {
    stop("all variables in the formula must be columns in data")
  }

  ## we have two possibilities
  ## 1. we have a left hand side, which is then the target
  ##    and the right hand side is used to define groups
  ##    --> perform multifactorial analysis of all pairs of groups
  ## 2. we have no left hand side
  ##    --> perform a mean analysis

  ## got a left hand side
  if (length(formula) == 3L){
    ## if an intercept is present remove it
    if (attr(terms(formula), "intercept") == 1)
      formula <- update(formula, . ~ . -1)
    ## get the model matrix
    modelMatrix <- model.matrix(formula, data = data)
    ## get the target
    yx <- eval(attr(terms(formula), "variables"), envir = data)
  ## no left hand side
  } else{

  }


}




