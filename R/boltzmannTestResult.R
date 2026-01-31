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
#'

new_boltzmannTestResult <- function(
    statistic,
    iDivergence,
    degreesOfFreedom,
    sampleSize,
    alternativeDistribution,
    alternativeMoments,
    hypothesisDistribution,
    hypothesisMoments,
    testedMoments,
    dataName
){
  res <- list(
    statistic = as.numeric(statistic),
    iDivergence = as.numeric(iDivergence),
    degreesOfFreedom = as.integer(degreesOfFreedom),
    sampleSize = as.integer(sampleSize),
    alternativeDistribution = alternativeDistribution,
    alternativeMoments = alternativeMoments,
    hypothesisDistribution = hypothesisDistribution,
    hypothesisMoments = hypothesisMoments,
    testedMoments = testedMoments,
    dataName = data.Name
  )
  class(res) = "boltzmannTestResult"
  validate_boltzmannTestResult(res)
}

validate_boltzmanntestResult(object){
  with(
    object,
    {
      if (! is.numeric(statistic) || length(statistic) != 1){
        stop("`statistic` must be a single numeric value")
      }
      if (! is.numeric(iDivergence) || length(iDivergence) != 1){
        stop("`iDivergence` must be a single numeric value")
      }
      if (!is.integer(degreesOfFreedom) || length(degreesOfFreedom) != 1){
        stop("`degreesOfFreedom` must be a single integer value")
      }
      if (dof < 1){
        stop("`degreesOfFreedom` must be at least 1")
      }
      if (!is.integer(sampleSize) || length(sampleSize) != 1){
        stop("`sampleSize` must be a single integer value")
      }
      if (sampleSize < 2){
        stop("`sampleSize` must be at least 2")
      }
      if (!is.numeric(alternativeDistribution) || !is.atomic(alternativeDistribution)){
        stop("`alternativeDistribution` must be a numeric vector")
      }
      if (!is.numeric(alternativeMoments) || !is.atomic(alternativeMoments)){
        stop("`alternativeMoments` must be a numeric vector")
      }
      if (!is.numeric(hypothesisDistribution) || !is.atomic(hypothesisDistribution)){
        stop("`hypothesisDistribution` must be a numeric vector")
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
    }
  )
  object
}
