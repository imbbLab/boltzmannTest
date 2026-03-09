#' Generic interface for the Boltzmann Test
#'
#' `boltzmann.test` with `object` being an outcome_tibble is the general
#' interface to perform the Boltzmann Test.
#'
#' @param outcomes an `outcomes_tibble` object
#' @param G a numeric matrix with the function values for the outcomes in the
#' rows and the outcomes in the columns for the hypothesis
#' @param eta a numeric vector with the same number of entries as `G` has rows
#' @param testedExpectations an integer vector indicating the tested
#' expectations
#' @param ambientExpectations an integer vector indicating the ambient
#' expectations
#' @param maxit an integer specifying the maximimal number of iterations used in
#' the I-projector (default 10000L)
#' @param tolerance a numeric specifying the numeric tolerance (default
#' .Machine$double.eps)
#'
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
#' @examples
#' ## nested hypothesis testing, deconfounding
#' ## generate outcome_tibble object
#' outcomes <- outcomes_tibble(kidneyStones)
#' ## coefficient matrix
#' G <- with(
#'   outcomes,
#'   rbind(
#'     norm = 1,
#'     # structural expectation fraction treatment A
#'     treatment_A = treatment == "A",
#'     # ambient expectation fraction small stone size given treatment B
#'     stoneSize_small.treatment_B =
#'       (stoneSize == "small" & treatment == "B") / 0.5,
#'     # ambient expectation fraction small stone size given treatment A
#'     stoneSize_small.treatment_A =
#'       (stoneSize == "small" & treatment == "A") / 0.5,
#'     # tested expectation: success rate treatment A versus B
#'     "success_yes:B_vs_A" = c(0, -2, 0, 2, 0, -2, 0, 2)
#'   )
#' )
#' ## hypothesized values of the expectations
#' eta <- c(
#'   norm = 1,
#'   treatment_A = 0.5,
#'   # overall frequency of small stones is 0.51
#'   stoneSize_small.treatment_B = 0.51,
#'   stoneSize_small.treatment_A = 0.51,
#'   "success_yes:B_vs_A" = 0.0
#' )
#' ## perform nested Boltzmann Test
#' (bt <- boltzmannTest(
#'   outcomes = outcomes,
#'   G = G,
#'   eta = eta,
#'   testedExpectations = 5,
#'   ambientExpectations = c(3,4)
#' ))
#' ## 1. project to a hypothesis with stone size frequencies the same in the
#' ##    treatment groups and not difference in the success rate between
#' ##    treatments
#' ## 2. project to the ambient alternative with stone size frequencies as the
#' ##    observed ones. Note the difference in the success rate in treatment B
#' ##    versus A of roughly 10.1%.
#' ##
#' ##    Sampling from a "population", where the treatments received the
#' ##    same number of patients with small and large stones with otherwise
#' ##    no overall difference in the treatment success rate, and keeping only
#' ##    samples with stone size frequencies equal to the observed ones
#' ##    exaggerates the difference between the treatments. Of course this is
#' ##    due to the dependency of success on the stone size.
#' boltzmann.test(success ~ stoneSize, data = kidneyStones)
#' ##    Because more patients with small stones received treatment B and
#' ##    because treating patients with small stones has a higher success rate,
#' ##    treatment B appears to have a higher success rate than treatment B even
#' ##   when there is no difference between the treatments.
#' ## 3. project back to the empirical. Now the difference in the success rate
#' ##    in treatment B versus A is 4.6%.
#' ## Since changing the frequencies of the stone sizes to the observed ones
#' ## already induces a shift in the difference from 0% to 10.1%, we see
#' ## that the real difference in the success rate is 4.6% - 10.1% = -5.5%.
#' ## The effect size estimate has reversed its sign as suggested by the
#' ## analysis stratified by stone size. Note, that this is not significant at
#' ## significance level 5%.
#'
#' @importFrom stats complete.cases
#' @importFrom stats pchisq
#' @importFrom tibble tibble
#' @export
#'
boltzmannTest <-function(
    outcomes, G, eta, testedExpectations, ambientExpectations = NULL,
    maxit = 10000L, tolerance = .Machine$double.eps){

  if (!is.matrix(G) || !is.numeric(G)){
    stop("`G` is not a numeric matrix")
  }
  if (!is.numeric(eta)){
    stop("`eta` must be numeric")
  }

  if (NROW(G) != length(eta)){
    stop("the number of generalized moments given by `eta` ",
         "does not match the number of rows in `G`")
  }
  if (any(!is.finite(eta))){
    stop("`eta` contains non-finite or missing values")
  }
  testedExpectations <- as.integer(testedExpectations)

  if (!is.null(ambientExpectations)){
    ambientExpectations <- as.integer(ambientExpectations)
    if (any(ambientExpectations < 1) || any(ambientExpectations > NROW(G))){
      stop("the indices of the ambient expectations are out of range")
    }
    if (any(ambientExpectations %in% testedExpectations)){
      stop("the ambient expectations must not include the tested expectations")
    }
  }


  if (any(testedExpectations < 1) || any(testedExpectations > NROW(G))){
    stop("the indices of the tested moments are out of range")
  }

  degreesOfFreedom <- length(testedExpectations)
  dataName <- deparse(substitute(outcomes))
  sampleSize <- sampleSize(outcomes)
  ## determine the sample expectations
  mu <- (G %*% empirical(outcomes))[, 1]

  etaNested <- NULL
  if (!is.null(ambientExpectations)){
    etaNested <- eta
    etaNested[ambientExpectations] <- mu[ambientExpectations]
  }

  ## 1. project empirical distribution to hypothesis linear family
  h <- iProjector(
    G = G,
    eta = eta,
    v = empirical(outcomes),
    maxit = maxit,
    convTolerance = tolerance
  )
  ## check for convergence
  if(h$converged != 0){
    ## here we cover the case, when some p's got zero
    if(any(h$p < tolerance)){
      warning("hypothesis conditions led to zero probabilities")
      return(boltzmannTestResult(
        statistic = Inf,
        iDivergence = Inf,
        degreesOfFreedom = degreesOfFreedom,
        sampleSize = sampleSize,
        pValue = 0,
        alternativeDistribution = empirical(outcomes),
        alternativeExpectations = mu,
        hypothesisDistribution = h$p,
        hypothesisExpectations = eta,
        testedExpectations = testedExpectations,
        dataName = dataName,
        coefficientMatrix = G,
        ambientDistribution =
          if (is.null(ambientExpectations)) NULL else rep(NA, NROW(outcomes)),
        ambientExpectations =
          if (is.null(ambientExpectations)) NULL else rep(NA, NROW(G))
      ))
    }
    ## In nested hypothesis testing we project first from the
    ## hypothesis distribution to the ambient alternative family
    baseDistribution <- if (!is.null(ambientExpectations)){
      p <- iProjector(
        G = G[-testedExpectations, ],
        eta = etaNested[-testedExpectations],
        v = h$p,
        maxit = maxit,
        convTolerance = tolerance
      )
      etaNested <- (G %*% p$p)[,1]
      p
    } else{
      h
    }
    iDivergence <- iDivergence(empirical(outcomes), baseDistribution$p)

    statistic <- 2 * sampleSize * iDivergence
    pValue <- pchisq(statistic, df = degreesOfFreedom, lower.tail = FALSE)
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
      ambientDistribution =
        if (is.null(ambientExpectations)) NULL else baseDistribution$p,
      ambientExpectations = etaNested
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
      ambientDistribution =
        if (is.null(ambientExpectations)) NULL else rep(NA, NROW(outcomes)),
      ambientExpectations =
        if (is.null(ambientExpectations)) NULL else etaNested
    )
  }




}
