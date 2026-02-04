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
#' Gives the number of tested generalized Moments and must be > 0
#' @param sampleSize a single integer value for the sample size. Must be > 1
#' @param pValue a single numeric value between 0 and 1 for the p-value
#' @param alternativeDistribution a numeric vector for the alternative distribution (usually equal
#' to the empirical distribution)
#' @param alternativeMoments a numeric vector for the values of the generalized moments
#' of the alternative distribution
#' @param hypothesisDistribution a numeric vector for the hypothesis distribution
#' @param hypothesisMoments a numeric vector for the values of the generalized moments
#' of the hypothesis distribution
#' @param testedMoments integer vector with the indices of the tested moments
#' @param dataName string with the name of the data
#' @param coefficientMatrix a numeric matrix holding the coefficient matrix to calculate
#' the generalized moments. The number of rows must equal the number of generalized moments and the number
#' of columns the number of entities.
#' @param nestedAlternativeDistribution a numeric vector for the nested alternative distribution
#' @param nestedAlternativeMoments a numeric vector for the values of the generalized moments
#' of the nested alternative distribution
#'
#' @export
new_boltzmannTestResult <- function(
    statistic,
    iDivergence,
    degreesOfFreedom,
    sampleSize,
    pValue,
    alternativeDistribution,
    alternativeMoments,
    hypothesisDistribution,
    hypothesisMoments,
    testedMoments,
    dataName,
    coefficientMatrix,
    method,
    nestedAlternativeDistribution = NULL,
    nestedAlternativeMoments = NULL
){
  if(!is.null(nestedAlternativeDistribution)){
    nestedAlternativeDistribution <- as.numeric(nestedAlternativeDistribution)
  }
  if(!is.null(nestedAlternativeMoments)){
    nestedAlternativeMoments <- as.numeric(nestedAlternativeMoments)
  }
  res <- list(
    statistic = statistic,
    iDivergence = iDivergence,
    degreesOfFreedom = degreesOfFreedom,
    sampleSize = sampleSize,
    pValue = pValue,
    alternativeDistribution = alternativeDistribution,
    alternativeMoments = alternativeMoments,
    hypothesisDistribution = hypothesisDistribution,
    hypothesisMoments = hypothesisMoments,
    testedMoments = testedMoments,
    dataName = dataName,
    coefficientMatrix = coefficientMatrix,
    method = method,
    nestedAlternativeDistribution = nestedAlternativeDistribution,
    nestedAlternativeMoments = nestedAlternativeDistribution
  )
  class(res) = "boltzmannTestResult"
  res
}
#' @description
#' `validate_boltzmannTestResult()` checks a `boltzmannTestResult` object for internal consistency.
#'
#' @param object a `boltzmannTestResult` class object
#' @rdname new_boltzmannTestResult
#' @export
validate_boltzmannTestResult <- function(object){
  ## check whether all list elements are present
  fields <- c(
    "statistic", "iDivergence", "degreesOfFreedom", "sampleSize",
    "pValue", "alternativeDistribution", "alternativeMoments",
    "hypothesisDistribution", "hypothesisMoments", "testedMoments",
    "dataName", "coefficientMatrix", "method", "nestedAlternativeDistribution",
    "nestedAlternativeMoments"
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
      if (!is.numeric(alternativeMoments) || !is.atomic(alternativeMoments)){
        stop("`alternativeMoments` must be a numeric vector")
      }
      if (!any(is.na(hypothesisDistribution))){
        if (!is.numeric(hypothesisDistribution) || !is.atomic(hypothesisDistribution)){
          stop("`hypothesisDistribution` must be a numeric vector")
        }
      }
      if (!is.numeric(hypothesisMoments) || !is.atomic(hypothesisMoments)){
        stop("`hypothesisMoments` must be a numeric vector")
      }
      if (length(hypothesisDistribution) != length(alternativeDistribution)){
        stop("`hypothesisDistribution` and `alternativeDistribution` must have the same length")
      }
      if (length(hypothesisMoments) != length(alternativeMoments)){
        stop("`hypothesisMoments` and `alternativeMoments` must have the same length")
      }
      if (!is.integer(testedMoments) ||!is.atomic(testedMoments)){
        stop("`testedMoments` must be a integer vector")
      }
      if (any(testedMoments < 1) || any(testedMoments > length(hypothesisMoments))){
        stop("`testedMoments` must be between 1 and the number of hypothesis moments")
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
      if (! is.null(nestedAlternativeMoments)){
        if (!is.numeric(nestedAlternativeMoments) || !is.atomic(nestedAlternativeMoments)){
          stop("`nestedAlternativeMoments` must be a numeric vector")
        }
      }

    }
  )
  object
}

#' @export
boltzmannTestResult <- function(
    statistic,
    iDivergence,
    degreesOfFreedom,
    sampleSize,
    pValue,
    alternativeDistribution,
    alternativeMoments,
    hypothesisDistribution,
    hypothesisMoments,
    testedMoments,
    dataName,
    coefficientMatrix,
    method = "Boltzmann Test",
    nestedAlternativeDistribution = NULL,
    nestedAlternativeMoments = NULL

){
  validate_boltzmannTestResult(
    new_boltzmannTestResult(
      statistic,
      iDivergence,
      degreesOfFreedom,
      sampleSize,
      pValue,
      alternativeDistribution,
      alternativeMoments,
      hypothesisDistribution,
      hypothesisMoments,
      testedMoments,
      dataName,
      coefficientMatrix,
      method,
      nestedAlternativeDistribution,
      nestedAlternativeMoments
    )
  )
}


#' @export
print.boltzmannTestResult <- function(object, digits = getOption("digits"), prefix = "\t", ...) {
  cat("\n")
  cat(strwrap(object$method, prefix = prefix), sep = "\n")
  cat("\n")
  cat("data:  ", object$dataName, "\n", sep = "")

  cat("statistic = ", object$statistic, ", df = ", object$degreesOfFreedom, sep = "")
  pval <- format.pval(object$pValue, digits = max(1L, digits - 3L))
  cat(", p-value ", if (substr(pval, 1L, 1L) == "<") pval else paste("=", pval), "\n", sep = "")

  testedMoments <- rep(FALSE, length(object$hypothesisMoments))
  testedMoments[object$testedMoments] = TRUE
  cat("hypothesis moments:\n")
  cat(
    ifelse(
      testedMoments,
      paste0("\033[0;", 1, "m",format(object$hypothesisMoments, digits = digits),"\033[0m"),
      paste0("\033[0;", 0, "m",format(object$hypothesisMoments, digits = digits),"\033[0m")
    )
  )
  cat("\n")
  if (!is.null(object$nestedAlternativeMoments)){
    cat("nested alternative moments:\n")
    cat(
      ifelse(
        testedMoments,
        paste0("\033[0;", 1, "m",format(object$nestedAlternativeMoments, digits = digits),"\033[0m"),
        paste0("\033[0;", 0, "m",format(object$nestedAlternativeMoments, digits = digits),"\033[0m")
      )
    )
    cat("\n")
  }
  cat("alternative (empirical) moments:\n")
  cat(
    ifelse(
      testedMoments,
      paste0("\033[0;", 1, "m",format(object$alternativeMoments, digits = digits),"\033[0m"),
      paste0("\033[0;", 0, "m",format(object$alternativeMoments, digits = digits),"\033[0m")
    )
  )
  cat("\ntested generalized moments are bold\n")
  cat("\n")

  invisible(object)
}
