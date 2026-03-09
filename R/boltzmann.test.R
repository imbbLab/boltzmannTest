#' Perform a Boltzmann Test for
#'
#' Performs a one-, two-, and multi-groups Boltzmann Test for the
#' the mean (one-group), the (paired) difference of the means (two-groups),
#' means and/or prevalences of multiple observables (class formula), means
#' and/or prevalences of multiple observables and multi-groups (class formula)
#'
#' @export
boltzmann.test <- function(x, ...) UseMethod("boltzmann.test")

#' @rdname boltzmann.test
#'
#' @param x a (non-empty) numeric vector of data values
#' @param y an optional (non-empty) numeric vector of data values
#' @param mu a single number indicating the true value of the mean or the
#' difference in means when performing a two sample test
#' @param paired a logical indicating whether to perform a paired test (default
#' FALSE)
#' @param ... further arguments passed to the different methods
#'
#' #' @returns
#' a list with class `boltzmannTestResult` containing the follow
#' components
#' \describe{
#'  \item{statistic}{the value of \eqn{\chi^2} statistic}
#'  \item{iDivergence}{the value of the I-divergence, a.k.a. Kullback-Leibler
#'  divergence of the hypothesis distribution or ambient alternative
#'  distribution to the empirical distribution}
#'  \item{degreesOfFreedom}{ the degrees of freedom used for computing the
#'  p-value via the \eqn{\chi^2} distribution with `degreesOfFreedom`}
#'  \item{sampleSize}{the sample size}
#'  \item{pValue}{the p-value}
#'  \item{alternativeDistribution}{the empirical distribution of outcomes}
#'  \item{alternativeExpectations}{the values of empirical expectations}
#'  \item{hypothesisDistribution}{the hypothesis distribution of outcomes}
#'  \item{hypothesisExpectations}{the values of hypothesized expectations}
#'  \item{testedExpectations}{an integer vector with the indices of the tested
#'    expectations}
#'  \item{dataName}{name of the input data}
#'  \item{coefficientMatrix}{an numeric matrix with the coefficient matrix}
#'  \item{method}{a string signifying the method}
#'  \item{ambientDistribution}{the ambient alternative distribution of outcomes
#'    if nested hypothesis testing is performed otherwise NULL}
#'  \item{ambientExpectations}{the values of the ambient expectations if nested
#'    hypothesis testing is performed otherwise NULL}
#' }
#'
#' @examples
#' ## examples using the NHANES data
#' data(nhanes)
#' ## testing whether females are on average the same height as males
#' with(
#'   nhanes,
#'   boltzmann.test(
#'     x = BMXHT[RIAGENDR == "female"],
#'     y = BMXHT[RIAGENDR == "male"]
#'    )
#' )
#'
#' ## testing whether females weigh on average 10 kg less than males
#' boltzmann.test(BMXWT ~ RIAGENDR, data = nhanes, nu = -10)
#'
#' ## testing whether females have the same average weight and height as males
#' boltzmann.test(cbind(BMXWT, BMXHT) ~ RIAGENDR, data = nhanes)
#'
#' data(kidneyStones)
#' ## testing whether the treatment success is the same in treatment A versus B
#' boltzmann.test(success ~ treatment, data = kidneyStones)
#'
#' ## testing whether the treatment success is the same in treatment A versus B
#' ## and small and large stone sizes (four groups)
#' boltzmann.test(success ~ treatment + stoneSize, data = kidneyStones)
#'
#' @importFrom stats complete.cases
#' @importFrom tibble tibble
#' @export
boltzmann.test.default <-function(x, y = NULL, mu = 0, paired = FALSE, ...){
  if (!missing(mu) && (length(mu) != 1 || is.na(mu)) || !is.numeric(mu)){
    stop("`mu` must be a single number")
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
        stop(
          paste0("`", xName, "` and `", yName, "` must have the same length")
        )
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
        group = factor(
          c(rep("x", sum(xok)), rep("y", sum(yok))),
          levels = c("x", "y")
        )
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
  res <- boltzmannTest(outcomes, G, eta, testedExpectations)
  res$dataName <- dataName
  res$method <- method
  res
}

#' @rdname boltzmann.test
#'
#' @param formula a formula either of the form `lhs` ~ `rhs`, where `lhs` is
#' either a numeric or factor variable or combinations thereof (see details) and
#' `rhs` is one or more factor variable used to define the groups, or of the
#' form ~ `rhs` where `rhs` contains one or more numeric or factor variable
#' (see details)
#' @param data a data.frame containing the variables used in `formula`
#' @param nu a single number or a numeric vector specifying the values of the
#' expectations
#' @details
#' The `formula` API allows to specify multi-group comparisons for one or more
#' target variable specified at `lhs`. If a target variable is numeric then the
#' difference of means of all groups except the first versus the mean of the
#' first group is compared. if a target variable is a factor then the difference
#' of the prevalence of the second (third, etc) level of that level of all
#' groups except the first versus the first group is compared.If only a `rhs`
#' is given, a multivariate comparison to the given `nu` is performed. If the
#' target is numeric the mean is calculated. If the target is a factor the
#' prevalence is calculated.
#'

#' @importFrom dplyr group_by across mutate cur_group_id ungroup
#' @importFrom stats model.frame complete.cases
#' @importFrom methods is
#' @export
boltzmann.test.formula <- function(formula, data, nu = 0, ...){

  ## coerce data to a data frame
  data <- as.data.frame(data)

  ## make sure that all variables used in the formula are also present
  ## as columns in the data
  vars <- all.vars(formula)
  if (!all(vars %in% colnames(data))) {
    stop("all variables in the formula must be columns in data")
  }

  varsClass <- vapply(
    vars,
    FUN = function(v) class(data[[v]]),
    FUN.VALUE = "character"
  )

  if (any(varsClass == "character")){
    stop("variables in `formula` are character vectors. ",
         "Convert these columns in `data` to factors")
  }
  completeCases <- completeCases_internal(data[, vars])
  if (any(completeCases == FALSE)){
    warning("selected columns in data contain unobserved values. ",
            "Removing the corresponding rows")
  }
  data <- subset(data, subset = completeCases, select = vars, drop = FALSE)

  ## we have two possibilities
  ## 1. we have a left hand side, which is then the target
  ##    and the right hand side is used to define groups
  ##    --> perform multifactorial analysis of all pairs of groups
  ## 2. we have no left hand side
  ##    --> perform a mean analysis

  ## got a left hand side

  if (length(formula) == 3L){

    ## get the model frame
    tryCatch(
      mf <- model.frame(
        formula = formula,
        data = data
      ),
      error = function(e) stop(paste0("formula is misspecified. ", e)),
      warning = function(e) stop(paste0("targets are misspecified. ", e))
    )

    ## get the target(s)
    y <- if (is.matrix(mf[[1]])){
      data.frame(mf[[1]])
    } else{
      mf[, 1, drop = FALSE]
    }

    for(targetVar in colnames(y)){
      if(targetVar %in% vars){
        if(is(data[[targetVar]], "factor")){
          y[[targetVar]] <- data[[targetVar]]
        }
      }
    }
    ## get the grouping variables
    X <- mf[,-1, drop = FALSE]

    xVarsClass <- vapply(
      colnames(X),
      FUN = function(v) class(data[[v]]),
      FUN.VALUE = "character"
    )

    if(any(xVarsClass != "factor")){
      stop("all grouping variables must be factors.")
    }
    ## make the groups

    groups <- group_by(X, across(everything()))
    groups <- mutate(groups, .group = cur_group_id())
    groupLabels <- lapply(
      attr(groups, "groups")[, colnames(X), drop = FALSE],
      as.character
    )
    groupLabels <- if(length(groupLabels) == 1){
      groupLabels[[1]]
    } else{
      groupLabels <- do.call(rbind, groupLabels)
      apply(
        groupLabels,
        2,
        paste,
        collapse = "."
      )
    }

    groups$.group <- factor(
      groups$.group,
      labels = groupLabels
    )
    groups <- ungroup(groups)
    groupLevels <- levels(groups$.group)
    if (length(groupLevels) < 2){
      stop("not enough groups to compare")
    }
    outcomes <- outcomes_tibble(
      data.frame(y, group = groups$.group, check.names = FALSE)
    )
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

    targetVars <- colnames(y)
  } else if (length(formula) == 2L){

    outcomes <- outcomes_tibble(data)
    G <- matrix(1, ncol = NROW(outcomes), dimnames = list("norm"))
    eta <- 1
    targetVars <- vars
    groupLevels <- NULL
    prevalences <- NULL
  }

  ## compute matrix elements for the expectations
  dt <- lapply(
    targetVars,
    function(targetVar){
      expectations(outcomes, targetVar, groupLevels, prevalences)
    }
  )
  ## target class
  targetVarClass <- vapply(
    targetVars,
    FUN = function(v) class(outcomes[[v]]),
    FUN.VALUE = "character"
  )
  ## the number of expectations per target variable
  numberExpectations <- vapply(dt, FUN = length, FUN.VALUE = 1)

  ## append to G
  dt <- do.call(rbind, unlist(dt, recursive = FALSE))

  targetExpectations <- seq_len(NROW(dt)) + NROW(G)
  G <- rbind(
    G,
    dt
  )
  targetExpectationsClass <- unlist(
    lapply(
      seq_along(targetVarClass),
      function(i) rep(targetVarClass[i], numberExpectations[i])
    )
  )

  testedExpectations <- if (length(nu) == 1){
    rep(nu, length(targetExpectations))
  } else{
    if (length(nu) != length(targetExpectations)){
      stop("the number of expectations in `nu` does not ",
           "match the number of tested moments")
    }
    nu
  }
  if (is.null(groupLevels)){
    if(any(
      targetExpectationsClass == "factor" &
      (testedExpectations <= 0 | testedExpectations >= 1))){
      stop("hypothesized prevalences outside (0, 1)")
    }
  }
  eta <- c(eta, testedExpectations)
  names(eta) <- rownames(G)

  boltzmannTest(outcomes, G, eta, targetExpectations)
}

#' @importFrom methods is
expectations <- function(
    outcomes, targetVar, groupLevels = NULL, prevalences = NULL){
  ## target is categorical
  if (is(outcomes[[targetVar]], "factor")){
    targetVarLevels <- levels(outcomes[[targetVar]])
    if (length(targetVarLevels) == 1){
      stop("categorical target ", targetVar, " has only one level")
    }
    res <- lapply(
      targetVarLevels[-1],
      function(targetVarLevel){
        if (is.null(groupLevels)){
          dt <- list(
            ifelse(
              outcomes[[targetVar]] == targetVarLevel,
              1,
              0
            )
          )
          names(dt) = paste0(
            targetVar, "_", targetVarLevel
          )
          dt
        } else{
          dt <-lapply(
            groupLevels[-1],
            function(group){
              ifelse(
                outcomes$group == group,
                (outcomes[[targetVar]] ==  targetVarLevel) / prevalences[group],
                ifelse(
                  outcomes$group == groupLevels[1],
                  -(outcomes[[targetVar]] == targetVarLevel) /
                    prevalences[groupLevels[1]],
                  0
                )
              )
            }
          )
          names(dt) <- paste0(
            targetVar, "_", targetVarLevel, ":",
            groupLevels[-1], "_vs_", groupLevels[1]
          )
          do.call(rbind, dt)
        }
      }
    )
    ## target is numeric
  } else{
    if (is.null(groupLevels)){
      dt <- list(outcomes[[targetVar]])
      names(dt) = targetVar
      dt
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
      names(dt) <- paste0(
        targetVar, ":", groupLevels[-1], "_vs_", groupLevels[1]
      )
      dt
    }

  }
}







