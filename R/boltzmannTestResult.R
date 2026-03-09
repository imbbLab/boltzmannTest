## boltzmannTestResult class

#' Constructor and Validator for the `boltzmannTestResult` class
#'
#' @description
#' Creates a `boltzmannTestResult` object
#'
#' @param statistic a single numeric value corresponding to the value of the
#' (test) statistic
#' @param iDivergence a single numeric value corresponding to the value of the
#' I-Divergence of the hypothesis distribution or ambient alternative
#' distribution to the (inner) alternative distribution
#' @param degreesOfFreedom a single integer value for the degrees of freedom
#' of the Chi^2 distribution. Gives the number of tested Expectations and
#' must be > 0
#' @param sampleSize a single integer value for the sample size. Must be > 1
#' @param pValue a single numeric value between 0 and 1 for the p-value
#' @param alternativeDistribution a numeric vector for the alternative
#' distribution (usually equal to the empirical distribution)
#' @param alternativeExpectations a numeric vector for the values of the
#' expectations of the alternative distribution
#' @param hypothesisDistribution a numeric vector for the
#' hypothesis distribution
#' @param hypothesisExpectations a numeric vector for the values of the
#' generalized moments of the hypothesis distribution
#' @param testedExpectations integer vector with the indices of the tested
#' expectations
#' @param dataName string with the name of the data
#' @param coefficientMatrix a numeric matrix holding the coefficient matrix
#' to calculate the generalized moments. The number of rows must equal the
#' number of generalized moments and the number of columns the number of
#' entities.
#' @param ambientDistribution a numeric vector for the ambient alternative
#' distribution
#' @param ambientExpectations a numeric vector for the values of the
#' expectations of the ambient alternative distribution
#'
#' @keywords internal
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
    ambientDistribution = NULL,
    ambientExpectations = NULL
){
  if(!is.null(ambientDistribution)){
    ambientDistribution <- as.numeric(ambientDistribution)
  }
  if(!is.null(ambientExpectations)){
    ambientExpectations <- as.numeric(ambientExpectations)
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
    ambientDistribution = ambientDistribution,
    ambientExpectations = ambientExpectations
  )
  class(res) <- "boltzmannTestResult"
  res
}
#' Validate a `boltzmannTestResult` object
#' @description
#' Checks a `boltzmannTestResult` object for internal consistency.
#'
#' @param object a `boltzmannTestResult` class object
#'
#' @keywords internal
validate_boltzmannTestResult <- function(object){
  ## check whether all list elements are present
  fields <- c(
    "statistic", "iDivergence", "degreesOfFreedom", "sampleSize",
    "pValue", "alternativeDistribution", "alternativeExpectations",
    "hypothesisDistribution", "hypothesisExpectations", "testedExpectations",
    "dataName", "coefficientMatrix", "method", "ambientDistribution",
    "ambientExpectations"
  )
  if (! all(fields %in% names(object))){
    stop("missing fields. These list elements are missing",
         paste(setdiff(fields, names(object)), collapse = ","))
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

      if (!is.numeric(alternativeDistribution)){
        stop("`alternativeDistribution` must be a numeric vector")
      }
      if (!is.numeric(alternativeExpectations)){
        stop("`alternativeExpectations` must be a numeric vector")
      }
      if (!any(is.na(hypothesisDistribution))){
        if (!is.numeric(hypothesisDistribution)){
          stop("`hypothesisDistribution` must be a numeric vector")
        }
      }
      if (!is.numeric(hypothesisExpectations)){
        stop("`hypothesisExpectations` must be a numeric vector")
      }
      if (length(hypothesisDistribution) != length(alternativeDistribution)){
        stop("`hypothesisDistribution` and `alternativeDistribution` ",
             "must have the same length")
      }
      if (length(hypothesisExpectations) != length(alternativeExpectations)){
        stop("`hypothesisExpectations` and `alternativeExpectations` ",
             "must have the same length")
      }
      if (!is.integer(testedExpectations) ||!is.atomic(testedExpectations)){
        stop("`testedExpectations` must be a integer vector")
      }
      if (any(
        testedExpectations < 1) ||
        any(testedExpectations > length(hypothesisExpectations))){
        stop("`testedExpectations` must be between 1 ",
          "and the number of hypothesis expectations")
      }

      if(!is.character(dataName)){
        stop("`dataName` must be a string")
      }

      if (!is.matrix(coefficientMatrix) || !is.numeric(coefficientMatrix)){
        stop("`coefficientMatrix` must be a numeric matrix")
      }
      if (! is.null(ambientDistribution)){
        if (!any(is.na(ambientDistribution))){
          if (!is.numeric(ambientDistribution)){
            stop("`ambientDistribution` must be a numeric vector")
          }
          if (length(hypothesisDistribution) != length(ambientDistribution)){
            stop("`hypothesisDistribution` and `ambientDistribution` ",
                 "must have the same length")
          }
        }

      }
      if (! is.null(ambientExpectations)){
        if (!is.numeric(ambientExpectations)){
          stop("`ambientExpectations` must be a numeric vector")
        }
        if (length(hypothesisExpectations) != length(ambientExpectations)){
          stop("`hypothesisExpectations` and `ambientExpectations` ",
               "must have the same length")
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
#' @param statistic a single numeric value corresponding to the value of the
#' (test) statistic
#' @param iDivergence a single numeric value corresponding to the value of the
#' I-Divergence of the hypothesis distribution or outer alternative
#' distribution to the (inner) alternative distribution
#' @param degreesOfFreedom a single integer value for the degrees of freedom
#' of the Chi^2 distribution. Gives the number of tested generalized
#' Expectations and must be > 0
#' @param sampleSize a single integer value for the sample size. Must be > 1
#' @param pValue a single numeric value between 0 and 1 for the p-value
#' @param alternativeDistribution a numeric vector for the alternative
#' distribution (usually equal to the empirical distribution)
#' @param alternativeExpectations a numeric vector for the values of the
#' observed expectations of the alternative distribution
#' @param hypothesisDistribution a numeric vector for the hypothesis
#' distribution
#' @param hypothesisExpectations a numeric vector for the values of the
#' expectations of the hypothesis distribution
#' @param testedExpectations integer vector with the indices of the
#' tested expectations
#' @param dataName string with the name of the data
#' @param coefficientMatrix a numeric matrix holding the coefficient matrix
#' to calculate the expectations. The number of rows must equal the number of
#' expectations and the number of columns the number of entities.
#' @param ambientDistribution a numeric vector for the ambient alternative
#' distribution
#' @param ambientExpectations a numeric vector for the values of the
#' expectations of the ambient alternative distribution
#' @returns
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
#'
#' @keywords internal
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
    ambientDistribution = NULL,
    ambientExpectations = NULL

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
      ambientDistribution,
      ambientExpectations
    )
  )
}


#' @export
print.boltzmannTestResult <- function(
    x, digits = getOption("digits"), prefix = "\t", ...) {
  object <- validate_boltzmannTestResult(x)
  cat("\n")
  cat(strwrap(object$method, prefix = prefix), sep = "\n")
  cat("\n")
  cat("data:  ", object$dataName, "\n", sep = "")

  cat(
    "statistic = ", object$statistic,
    ", df = ", object$degreesOfFreedom,
    ", N = ", object$sampleSize, sep = ""
  )
  pval <- format.pval(object$pValue, digits = max(1L, digits - 3L))
  cat(", p-value ", if (substr(pval, 1L, 1L) == "<") pval else
    paste("=", pval), "\n\n", sep = "")

  testedExpectations <- rep(FALSE, length(object$hypothesisExpectations))
  testedExpectations[object$testedExpectations] <- TRUE
  moments <- data.frame(
    hypothesis = object$hypothesisExpectations,
    alternative = object$alternativeExpectations
  )

  if (!is.null(object$ambientExpectations)){
    moments$ambient <- object$ambientExpectations
    moments <- moments[, c("hypothesis", "ambient", "alternative")]
  }
  moments$tested <- ifelse(testedExpectations, "*", "")
  print(moments, digits = digits)

  invisible(object)
}
