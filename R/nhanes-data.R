#' NHANES 2021-2023 data
#'
#' Data from NHANES cycle 2021-2023
#'
#' @docType data
#'
#' @usage data(kidneyStones)
#'
#' @format A data.frame with 5620 observations of five variables (SEQN,
#' RIAGENDR, WTMEC2YR, BMXWT, BMXHT)
#'
#' @keywords dataset
#'
#' @references
#' https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2021/DataFiles/DEMO_L.xpt
#' https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2021/DataFiles/BMX_L.xpt
#'
#' @examples
#' data(nhanes)
#' ## testing whether the average BMI is equal in females vs males
#' ## no reweighting by WTMEC2YR
#' boltzmann.test(BMXWT/(BMXHT / 100)^2 ~ RIAGENDR, data = nhanes)
"nhanes"
