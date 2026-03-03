## boltzmannTestResult class

#' Constructor and Validator for the `boltzmannTestResult` class
#'
#' @description
#' Creates a `boltzmannTestResult` object
#'
#' @param statistic a single numeric value corresponding to the value of the (test) statistic
#' @param iDivergence a single numeric value corresponding to the value of the I-Divergence
#' of the hypothesis distribution or outer alternative distribution to the (inner) alternative
#' distribution
#' @param degreesOfFreedom a single integer value for the degrees of freedom of the Chi^2 distribution.
#' Gives the number of tested generalized Expectations and must be > 0
#' @param sampleSize a single integer value for the sample size. Must be > 1
#' @param pValue a single numeric value between 0 and 1 for the p-value
#' @param alternativeDistribution a numeric vector for the alternative distribution (usually equal
#' to the empirical distribution)
#' @param alternativeExpectations a numeric vector for the values of the generalized moments
#' of the alternative distribution
#' @param hypothesisDistribution a numeric vector for the hypothesis distribution
#' @param hypothesisExpectations a numeric vector for the values of the generalized moments
#' of the hypothesis distribution
#' @param testedExpectations integer vector with the indices of the tested moments
#' @param dataName string with the name of the data
#' @param coefficientMatrix a numeric matrix holding the coefficient matrix to calculate
#' the generalized moments. The number of rows must equal the number of generalized moments and the number
#' of columns the number of entities.
#' @param nestedAlternativeDistribution a numeric vector for the nested alternative distribution
#' @param nestedAlternativeExpectations a numeric vector for the values of the generalized moments
#' of the nested alternative distribution
#' @noRd
#' @export
new_boltzmannTestResult <- function(
    statistic,
    iDivergence,
    degreesOfFreedom,
    sampleSize,
    pValue,
    alternativeDistribution,
    alternativeExpectations,
    hypothesisDistribution,
    hypothesisExpectations,
    testedExpectations,
    dataName,
    coefficientMatrix,
    method,
    nestedAlternativeDistribution = NULL,
    nestedAlternativeExpectations = NULL
){
  if(!is.null(nestedAlternativeDistribution)){
    nestedAlternativeDistribution <- as.numeric(nestedAlternativeDistribution)
  }
  if(!is.null(nestedAlternativeExpectations)){
    nestedAlternativeExpectations <- as.numeric(nestedAlternativeExpectations)
  }
  res <- list(
    statistic = statistic,
    iDivergence = iDivergence,
    degreesOfFreedom = degreesOfFreedom,
    sampleSize = sampleSize,
    pValue = pValue,
    alternativeDistribution = alternativeDistribution,
    alternativeExpectations = alternativeExpectations,
    hypothesisDistribution = hypothesisDistribution,
    hypothesisExpectations = hypothesisExpectations,
    testedExpectations = testedExpectations,
    dataName = dataName,
    coefficientMatrix = coefficientMatrix,
    method = method,
    nestedAlternativeDistribution = nestedAlternativeDistribution,
    nestedAlternativeExpectations = nestedAlternativeDistribution
  )
  class(res) = "boltzmannTestResult"
  res
}
#' @description
#' `validate_boltzmannTestResult()` checks a `boltzmannTestResult` object for internal consistency.
#'
#' @param object a `boltzmannTestResult` class object
#' @noRd
#' @export
validate_boltzmannTestResult <- function(object){
  ## check whether all list elements are present
  fields <- c(
    "statistic", "iDivergence", "degreesOfFreedom", "sampleSize",
    "pValue", "alternativeDistribution", "alternativeExpectations",
    "hypothesisDistribution", "hypothesisExpectations", "testedExpectations",
    "dataName", "coefficientMatrix", "method", "nestedAlternativeDistribution",
    "nestedAlternativeExpectations"
  )
  if (! all(fields %in% names(object))){
    stop("these list elements are missing", setdiff(fields, names(object)))
  }

  with(
    object,
    {
      if (!any(is.na(statistic))){
        if (!is.numeric(statistic) || length(statistic) != 1){
          stop("`statistic` must be a single numeric value")
        }
      }
      if (!any(is.na(iDivergence))){
        if (! is.numeric(iDivergence) || length(iDivergence) != 1){
          stop("`iDivergence` must be a single numeric value")
        }
      }
      if (!any(is.na(degreesOfFreedom))){
        if (!is.integer(degreesOfFreedom) || length(degreesOfFreedom) != 1){
          stop("`degreesOfFreedom` must be a single integer value")
        }
      }
      if (degreesOfFreedom < 1){
        stop("`degreesOfFreedom` must be at least 1")
      }
      if (!is.integer(sampleSize) || length(sampleSize) != 1){
        stop("`sampleSize` must be a single integer value")
      }
      if (sampleSize < 2){
        stop("`sampleSize` must be at least 2")
      }
      if (!any(is.na(pValue))){
        if (!is.numeric(pValue) || length(pValue) != 1){
          stop("`pValue` must be a single numeric value")
        }
        if (pValue < 0 || pValue > 1){
          stop("`pValue` must be between 0 and 1")
        }
      }

      if (!is.numeric(alternativeDistribution) || !is.atomic(alternativeDistribution)){
        stop("`alternativeDistribution` must be a numeric vector")
      }
      if (!is.numeric(alternativeExpectations) || !is.atomic(alternativeExpectations)){
        stop("`alternativeExpectations` must be a numeric vector")
      }
      if (!any(is.na(hypothesisDistribution))){
        if (!is.numeric(hypothesisDistribution) || !is.atomic(hypothesisDistribution)){
          stop("`hypothesisDistribution` must be a numeric vector")
        }
      }
      if (!is.numeric(hypothesisExpectations) || !is.atomic(hypothesisExpectations)){
        stop("`hypothesisExpectations` must be a numeric vector")
      }
      if (length(hypothesisDistribution) != length(alternativeDistribution)){
        stop("`hypothesisDistribution` and `alternativeDistribution` must have the same length")
      }
      if (length(hypothesisExpectations) != length(alternativeExpectations)){
        stop("`hypothesisExpectations` and `alternativeExpectations` must have the same length")
      }
      if (!is.integer(testedExpectations) ||!is.atomic(testedExpectations)){
        stop("`testedExpectations` must be a integer vector")
      }
      if (any(testedExpectations < 1) || any(testedExpectations > length(hypothesisExpectations))){
        stop("`testedExpectations` must be between 1 and the number of hypothesis moments")
      }

      if(!is.character(dataName)){
        stop("`dataName` must be a string")
      }

      if (!is.matrix(coefficientMatrix) || !is.atomic(coefficientMatrix)){
        stop("`coefficientMatrix` must be a numeric matrix")
      }
      if (! is.null(nestedAlternativeDistribution)){
        if (!any(is.na(nestedAlternativeDistribution))){
          if (!is.numeric(nestedAlternativeDistribution) || !is.atomic(nestedAlternativeDistribution)){
            stop("`nestedAlternativeDistribution` must be a numeric vector")
          }
        }
      }
      if (! is.null(nestedAlternativeExpectations)){
        if (!is.numeric(nestedAlternativeExpectations) || !is.atomic(nestedAlternativeExpectations)){
          stop("`nestedAlternativeExpectations` must be a numeric vector")
        }
      }

    }
  )
  object
}

#' Build a `boltzmannTestResult` class object
#'
#' @description
#' Creates a `boltzmannTestResult` object
#'
#' @param statistic a single numeric value corresponding to the value of the (test) statistic
#' @param iDivergence a single numeric value corresponding to the value of the I-Divergence
#' of the hypothesis distribution or outer alternative distribution to the (inner) alternative
#' distribution
#' @param degreesOfFreedom a single integer value for the degrees of freedom of the Chi^2 distribution.
#' Gives the number of tested generalized Expectations and must be > 0
#' @param sampleSize a single integer value for the sample size. Must be > 1
#' @param pValue a single numeric value between 0 and 1 for the p-value
#' @param alternativeDistribution a numeric vector for the alternative distribution (usually equal
#' to the empirical distribution)
#' @param alternativeExpectations a numeric vector for the values of the generalized moments
#' of the alternative distribution
#' @param hypothesisDistribution a numeric vector for the hypothesis distribution
#' @param hypothesisExpectations a numeric vector for the values of the generalized moments
#' of the hypothesis distribution
#' @param testedExpectations integer vector with the indices of the tested moments
#' @param dataName string with the name of the data
#' @param coefficientMatrix a numeric matrix holding the coefficient matrix to calculate
#' the generalized moments. The number of rows must equal the number of generalized moments and the number
#' of columns the number of entities.
#' @param nestedAlternativeDistribution a numeric vector for the nested alternative distribution
#' @param nestedAlternativeExpectations a numeric vector for the values of the generalized moments
#' of the nested alternative distribution
#' @export
boltzmannTestResult <- function(
    statistic,
    iDivergence,
    degreesOfFreedom,
    sampleSize,
    pValue,
    alternativeDistribution,
    alternativeExpectations,
    hypothesisDistribution,
    hypothesisExpectations,
    testedExpectations,
    dataName,
    coefficientMatrix,
    method = "Boltzmann Test",
    nestedAlternativeDistribution = NULL,
    nestedAlternativeExpectations = NULL

){
  validate_boltzmannTestResult(
    new_boltzmannTestResult(
      statistic,
      iDivergence,
      degreesOfFreedom,
      sampleSize,
      pValue,
      alternativeDistribution,
      alternativeExpectations,
      hypothesisDistribution,
      hypothesisExpectations,
      testedExpectations,
      dataName,
      coefficientMatrix,
      method,
      nestedAlternativeDistribution,
      nestedAlternativeExpectations
    )
  )
}


#' @export
print.boltzmannTestResult <- function(object, digits = getOption("digits"), prefix = "\t", ...) {
  cat("\n")
  cat(strwrap(object$method, prefix = prefix), sep = "\n")
  cat("\n")
  cat("data:  ", object$dataName, "\n", sep = "")

  cat("statistic = ", object$statistic, ", df = ", object$degreesOfFreedom, ", N = ", object$sampleSize, sep = "")
  pval <- format.pval(object$pValue, digits = max(1L, digits - 3L))
  cat(", p-value ", if (substr(pval, 1L, 1L) == "<") pval else paste("=", pval), "\n\n", sep = "")

  testedExpectations <- rep(FALSE, length(object$hypothesisExpectations))
  testedExpectations[object$testedExpectations] = TRUE
  moments <- data.frame(
    hypothesis = object$hypothesisExpectations,
    alternative = object$alternativeExpectations
  )

  if (!is.null(object$nestedAlternativeExpectations)){
    moments$nested <- object$nestedAlternativeExpectations
    moments <- moments[, c("hypothesis", "nested", "alternative")]
  }
  moments$tested <- ifelse(testedExpectations, "*", "")
  print(moments, digits = digits)

  invisible(object)
}
