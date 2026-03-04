#' Kidney stone data
#'
#' Data from Charig et al. (1986), Table 1. Treatment A combines Nephrolithotomy/pyelolithotomy,
#' Pyelolithotomy, Ureterolithotomy, treatment B corresponds to Percutaneous nephrolithotomy
#'
#' @docType data
#'
#' @usage data(kidneyStones)
#'
#' @format A data.frame with 700 observations of three factor variables (stoneSize,
#' treatment, success)
#'
#' @keywords dataset
#'
#' @references Charig et al. (1986) Comparison of treatment of renal
#' calculi by open surgery, percutaneous nephrolithotomy,
#' and extracorporeal shockwave lithotripsy. BMJ 292 (6524), 879-882
#'
#' @examples
#' data(kidneyStones)
#' ## testing treatment success of treatment B versus treatment A
#' boltzmann.test(success ~ treatment, data = kidneyStones)
#' ## testing treatment success for small versus large stones
#' boltzmann.test(success ~ stoneSize, data = kidneyStones)
"kidneyStones"
