#' Perform a Boltzmann Test
#' @export
boltzmann.test <- function(object, ...) UseMethod("boltzmann.test")

#' @rdname boltzmann.test
#' @description
#' `boltzmann.test` with `object` being an outcome_tibble is the general
#' interface to perform the Boltzmann Test.
#'
#' @param outcomes an `outcomes_tibble` object
#' @param G a numeric matrix with the function values for the outcomes in the rows
#' and the outcomes in the columns
#' @param eta a numeric vector with the same number of entries as `G` has rows
#' @param testedExpectations an integer vector indicating the tested expectations
#' @param nestedExpectations an integer vector indicating the nested expectations
#' @param maxit an integer specifying the maximimal number of iterations used in
#' the I-projector (default 10000L)
#' @param tolerance a numeric specifying the numeric tolerance (default
#' .Machine$double.eps)
#' @importFrom tibble tibble
#' @export
#'
boltzmann.test.outcomes_tibble <-function(outcomes, G, eta, testedExpectations, nestedExpectations = NULL, maxit = 10000L, tolerance = .Machine$double.eps){
  if (NROW(G) != length(eta)){
    stop("the number of generalized moments given by `eta` does not match the number of rows in `G`")
  }
  if (!is.null(nestedExpectations)){
    nestedExpectations <- as.vector(nestedExpectations)
    if (any(nestedExpectations) < 1 || any(nestedExpectations) > NROW(G)){
      stop("the indices of the nested expectations are out of range")
    }
    if (any(nestedExpectations %in% testedExpectations)){
      stop("the nested expectations must not include the tested expectations")
    }
  }
  testedExpectations <- as.integer(testedExpectations)
  if (any(testedExpectations < 1) || any(testedExpectations > NROW(G))){
    stop("the indices of the tested moments are out of range")
  }
  if (length(testedExpectations) < 1){
    stop("there must be at least one moment to be tested")
  }
  degreesOfFreedom <- length(testedExpectations)
  dataName <- deparse(substitute(outcomes))
  sampleSize <- sampleSize(outcomes)
  ## determine the generalized moments for the empirical distribution
  mu <- (G %*% empirical(outcomes))[, 1]

  etaNested <- NULL
  if (!is.null(nestedExpectations)){
    etaNested <- eta
    eta[nestedExpectations] = mu[nestedExpectations]
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
        alternativeExpectations = mu,
        hypothesisDistribution = rep(NA, NROW(outcomes)),
        hypothesisExpectations = eta,
        testedExpectations = testedExpectations,
        dataName = dataName,
        coefficientMatrix = G,
        nestedAlternativeDistribution = if (is.null(nestedExpectations)) NULL else rep(NA, NROW(outcomes)),
        nestedAlternativeExpectations = if (is.null(nestedExpectations)) NULL else etaNested
      ))
      ## we may remove the offending entry and try again?
    }
    ## "update" hypothesis distribution by the empirical moments

    baseDistribution <- if (!is.null(nestedExpectations)){
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
      alternativeExpectations = mu,
      hypothesisDistribution = h$p,
      hypothesisExpectations = eta,
      testedExpectations = testedExpectations,
      dataName = dataName,
      coefficientMatrix = G,
      nestedAlternativeDistribution = if (is.null(nestedExpectations)) NULL else baseDistribution$p,
      nestedAlternativeExpectations = etaNested
    )
  } else{
    boltzmannTestResult(
      statistic = NA,
      iDivergence = NA,
      degreesOfFreedom = degreesOfFreedom,
      sampleSize = sampleSize,
      pValue = NA,
      alternativeDistribution = empirical(outcomes),
      alternativeExpectations = mu,
      hypothesisDistribution = rep(NA, NROW(outcomes)),
      hypothesisExpectations = eta,
      testedExpectations = testedExpectations,
      dataName = dataName,
      coefficientMatrix = G,
      nestedAlternativeDistribution = if (is.null(nestedExpectations)) NULL else rep(NA, NROW(outcomes)),
      nestedAlternativeExpectations = if (is.null(nestedExpectations)) NULL else etaNested
    )
  }




}
#' @rdname boltzmann.test
#' @description
#' `boltzmann.test` with `object` being a numeric vector performs one or two sample
#' Boltzmann Tests for the mean. It follows the `t.test` interface.
#'
#' @param x a (non-empty) numeric vector of data values
#' @param y an optional (non-empty) numeric vector of data values
#' @param mu a single number indicating the true value of the mean or the difference in
#' means when performing a two sample test
#' @param paired a logical indicating whether to perform a paired test (default
#' FALSE)
#'
#' @importFrom tibble tibble
#' @export
boltzmann.test.numeric <-function(x, y = NULL, mu = 0, paired = FALSE){
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
    dataName <- paste0(xName, " and ", yName)
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
      names(eta) <- c("norm", paste0(xName, "-", yName))
      rownames(G) <- names(eta)
      testedExpectations <- 2
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
          group == "x"
        ),
        with(
          outcomes,
          ifelse(
            group == "x",
            z / xNameProb,
            -z / (1 - xNameProb)
          )
        )

      )
      eta <- c(1, xNameProb, mu)
      names(eta) <- c("norm", xName, paste0(xName, "_vs_", yName))
      testedExpectations <- 3
      rownames(G) <- names(eta)
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
    names(eta) <- c("norm", xName)
    rownames(G) <- names(eta)
    testedExpectations <- 2
  }
  res <- boltzmann.test(outcomes, G, eta, testedExpectations)
  res$dataName = dataName
  res$method = method
  res
}

#' @rdname boltzmann.test
#' @description
#' `boltzmann.test` with `object` being a formula performs multigroup
#' or multivariate comparisons
#'
#' @param formula a formula either of the form `lhs` ~ `rhs`, where `lhs` is either
#' a numeric or factor variable or combinations thereof (see details) and
#' `rhs` is one or more factor variable used to define the groups, or of the form
#' ~ `rhs` where `rhs` contains one or more numeric or factor variable (see details)
#' @param data a data.frame containing the variables used in `formula`
#' @param nu a single number or a numeric vector specifying the values of the
#' expectations
#' @details
#' The `formula` API allows to specify multigroup comparisons for one or more target
#' variable specified at `lhs`. If a target variable is numeric then the difference
#' of means of all groups except the first versus the mean of the first group is compared.
#' if a target variable is a factor then the difference of the prevalence of the second (third, etc) level
#' of that level of all groups except the first versus the first group is compared.
#' If only a `rhs` is given, a multivariate comparison to the given `nu` is performed. If the target
#' is numeric the mean is calculated. If the target is a factor the prevalence is calculated.
#'

#' @importFrom dplyr group_by across mutate cur_group_id ungroup
#' @importFrom stats model.frame
#' @export
boltzmann.test.formula <- function(formula, data, nu = 0){
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

  varsClass <- sapply(vars, function(v) class(data[[v]]))

  if (any(varsClass == "character")){
    stop("variables in `formula` are character vectors. Convert these columns in `data` to factors")
  }

  if (any(complete.cases(data[vars]) == FALSE)){
    warning("selected columns in data contain unobserved values. Removing the corresponding rows")
    data = data[complete.cases(data[vars]), vars, drop = FALSE]
  } else{
    data = data[, vars, drop = FALSE]
  }

  ## we have two possibilities
  ## 1. we have a left hand side, which is then the target
  ##    and the right hand side is used to define groups
  ##    --> perform multifactorial analysis of all pairs of groups
  ## 2. we have no left hand side
  ##    --> perform a mean analysis

  ## got a left hand side

  if (length(formula) == 3L){


    tryCatch(
      mf <- model.frame(
        formula = formula,
        data = data
      ),
      error = function(e) stop(paste0("formula is misspecified. ", e)),
      warning = function(e) stop(paste0("`targets are misspecified. ", e))
    )

    ## get the target(s)
    y <- mf[, 1, drop = FALSE]
    if (is.matrix(y[[1]])){
      y <- data.frame(y[[1]])
    }
    for(targetVar in colnames(y)){
      if(targetVar %in% vars){
        if(class(data[[targetVar]]) == "factor"){
          y[[targetVar]] <- data[[targetVar]]
        }
      }
    }
    ## get the grouping variables
    X <- mf[,-1, drop = FALSE]

    xVarsClass <- sapply(colnames(X), function(v) class(X[[v]]))

    if(any(xVarsClass != "factor")){
      stop("all grouping variables must be factors")
    }
    ## make the groups

    groups <- group_by(X, across(everything()))
    groups <- mutate(groups, .group = cur_group_id())
    groups$.group = factor(
      groups$.group,
      labels = apply(
        sapply(
          attr(groups, "groups")[, colnames(X)],
          as.character
        ),
        1,
        paste,
        collapse = "."
      )
    )
    groups <- ungroup(groups)
    groupLevels <- levels(groups$.group)
    if (length(groupLevels) < 2){
      stop("not enough groups to compare")
    }
    outcomes <- outcomes_tibble(data.frame(y, group = groups$.group, check.names = FALSE))
    ## group prevalences except first
    prevalence <- lapply(
      groupLevels,
      function(group){
        ifelse(outcomes$group == group, 1, 0 )
      }
    )
    names(prevalence) <- groupLevels
    prevalence <- do.call(rbind, prevalence)

    G <- rbind(norm = 1, prevalence[-1, , drop = FALSE])

    prevalences <- (prevalence %*% empirical(outcomes))[, 1]
    eta <- c(norm = 1, prevalences[-1])

    ## all groups versus first
    for (targetVar in colnames(y)){
      if(class(outcomes[[targetVar]]) == "factor"){
        targetVarLevels <- levels(outcomes[[targetVar]])
        if (length(targetVarLevels) == 1){
          stop("categorical target ", targetVar, "has only one level")
        }
        for (targetVarLevel in targetVarLevels[-1]){
          dt <-lapply(
            groupLevels[-1],
            function(group){
              ifelse(
                outcomes$group == group,
                (outcomes[[targetVar]] ==  targetVarLevel) / prevalences[group],
                ifelse(
                  outcomes$group == groupLevels[1],
                  -(outcomes[[targetVar]] == targetVarLevel) / prevalences[groupLevels[1]],
                  0
                )

              )
            }
          )
          names(dt) <- paste0(targetVar, "_", targetVarLevel, ":", groupLevels[-1], "_vs_", groupLevels[1])
          dt <- do.call(rbind, dt)
          G <- rbind(G, dt)
        }

      } else{
        dt <-lapply(
          groupLevels[-1],
          function(group){
            ifelse(
              outcomes$group == group,
              outcomes[[targetVar]] / prevalences[group],
              ifelse(
                outcomes$group == groupLevels[1],
                -outcomes[[targetVar]] / prevalences[groupLevels[1]],
                0
              )

            )
          }
        )
        names(dt) <- paste0(targetVar, ":", groupLevels[-1], "_vs_", groupLevels[1])
        dt <- do.call(rbind, dt)
        G <- rbind(G, dt)
      }
    }
    targetExpectations <- seq(length(eta) + 1, nrow(G))
    if (length(nu) == 1){
      eta <- c(eta, rep(nu, length(targetExpectations)))
    } else{
      if (length(nu) != length(targetExpectations)){
        stop("the number of expectations in `nu` does not match the number of tested moments")
      }
      eta <- c(eta, nu)
    }
    names(eta) = rownames(G)

  } else if (length(formula) == 2){

    outcomes <- outcomes_tibble(data)
    nam <- NULL
    dt <- lapply(
      vars,
      function(targetVar){
        if (class(outcomes[[targetVar]]) == "factor"){
          targetVarLevels <- levels(outcomes[[targetVar]])
          if (length(targetVarLevels) == 1){
            stop("categorical target ", targetVar, "has only one level")
          }
          res <- lapply(
            targetVarLevels[-1],
            function(targetVarLevel){
              ifelse(
                outcomes[[targetVar]] == targetVarLevel,
                1,
                0
              )
            }
          )
          names(res) <- paste0(targetVar, "_", targetVarLevels[-1])
          do.call(rbind, res)
        } else{
          outcomes[[targetVar]]
        }
      }
    )

    targetVarClass <- sapply(vars, function(v) class(outcomes[[v]]))
    names(dt)[targetVarClass != "factor"] <- vars[targetVarClass != "factor"]
    numberExpectations <- sapply(dt, NROW)
    numberExpectations[targetVarClass != "factor"] <- 1
    checkClass <- unlist(lapply(seq_along(targetVarClass), function(i) rep(targetVarClass[i], numberExpectations[i])))


    dt <- do.call(rbind, dt)
    G <- rbind(norm = 1, dt)
    targetExpectations <- seq(2, nrow(G))
    if (length(nu) == 1){
      eta <- rep(nu, length(targetExpectations))
    } else{
      if (length(nu) != length(targetExpectations)){
        stop("the number of expectations in `nu` does not match the number of tested moments")
      }
      eta <- nu
    }
    if(any(checkClass == "factor" & (eta <= 0 | eta >= 1))){
      stop("hypothesized prevalences outside (0, 1)")
    }
    eta <- c(1, eta)
    names(eta) = rownames(G)
  } else{
    stop("unrecognized `formula`")
  }
  boltzmann.test(outcomes, G, eta, targetExpectations)
}




