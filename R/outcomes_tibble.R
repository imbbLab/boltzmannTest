

#' @importFrom dplyr group_by across everything summarise n select
#' @importFrom tibble new_tibble validate_tibble
#'
#' @noRd
new_outcomes_tibble <- function(data, empirical, sampleSize){
  if(! is.data.frame(data)){
    stop("`data` must be a data frame or tibble")
  }
  if(!is.numeric(empirical)){
    stop("`empirical` must be numeric")
  }
  if(length(empirical) != nrow(data)){
    stop("`empirical` must have one value per row of `data`. There are", length(empirical), "entries in `empirical` and", nrow(data), "rows in `data`")
  }
  if(length(sampleSize) != 1 || !is.numeric(sampleSize)){
    stop("`sampleSize` must be a single numeric value")
  }
  validate_outcomes_tibble(
    new_tibble(
      data,
      nrow = nrow(data),
      empirical = empirical,
      sampleSize = sampleSize,
      class = "outcomes_tibble"
    )
  )
}
#' @description
#' `validate_outcomes_tibble()` checks an `outcomes_tibble` for internal consistency.
#'
#' @param object an outcomes_tibble class object
#' @noRd
validate_outcomes_tibble <- function(object){
  object <- validate_tibble(object)
  empirical <- attr(object, "empirical", exact = TRUE)
  sampleSize <- attr(object, "sampleSize", exact = TRUE)
  if(is.null(empirical)){
    stop("`empirical` attribute is missing")
  }
  if(is.null(sampleSize)){
    stop("`sampleSize` attribute is missing")
  }
  if(!is.numeric(empirical)){
    stop("`empirical` must be numeric")
  }
  if(length(empirical) != nrow(object)){
    stop("`empirical` must have one value per row of `data`. There are", length(empirical), "entries in `empirical` and", nrow(object), "rows in `data`")
  }
  if (any(empirical < 0)){
    stop("`empirical` contains negative entries")
  }
  if (abs(sum(empirical) - 1) > sqrt(.Machine$double.eps)){
    stop("`empirical` does not sum to 1")
  }
  if (sampleSize < nrow(object)){
    stop("fewer observations", sampleSize, "than outcomes", nrow(object), call. = FALSE)
  }
  if (!is.numeric(sampleSize) || length(sampleSize) != 1) {
    stop("`sampleSize` must be a single number.")
  }
  if (sampleSize <= 0) {
    stop("`sampleSize` must be positive.")
  }
  object
}
#' Build an `outcomes_tibble` class object
#'
#' @description
#' Constructs an `outcomes_tibble()`  a subclass of `tibble` with unique rows from an input tibble-like object. It adds two attributes
#' \itemize{
#'   \item{`empirical` the empirical relative frequencies of the outcomes (unique rows)}
#'   \item{`sampleSize` the sample size given by the number of rows of the input `data.frame`}
#' }
#'
#' @param data a data frame
#' @export
outcomes_tibble <- function(data){
  if (!is.data.frame(data)){
    stop("df should be a data frame but is", class(data))
  }
  ## the number of rows is the sample size
  sampleSize <- nrow(data)
  ## group data by all columns
  data <- group_by(data, across(everything()))
  ## summarize and drop
  data <- summarise(data, empirical = n() / sampleSize, .groups = "drop")

  new_outcomes_tibble(
    data = select(data, -empirical),
    empirical = data$empirical,
    sampleSize = sampleSize
  )

}

#' `dplyr` row slicing
#'
#' @description
#' This function is automatically called after the `dplyr` functions `filter`, `slice`,
#' `arrange` etc. It automatically updates (marginalize) the empirical relative frequencies
#' and the sample size
#'
#' @param object an `outcomes_tibble` object
#' @param i row indices selected by an upstream `dplyr` function
#' @param ... additional parameters
#'
#' @importFrom dplyr dplyr_row_slice
#' @noRd
#' @export

dplyr_row_slice.outcomes_tibble <- function(object, i, ...) {
  if (length(i) == 0) {
    stop("subsetting removed all rows.")
  }

  out <- NextMethod()

  empirical <- attr(object, "empirical")[i]
  sampleSize  <- attr(object, "sampleSize")


  validate_outcomes(
    new_outcomes(
      out,
      empirical = empirical / sum(empirical),
      sampleSize = sum(empirical * sampleSize)
    )
  )
}

#' Subset an `outcomes_tibble`
#'
#' Subsetting an `outcomes_tibble` works similarily to subsetting a `tibble`,
#' but the object maintains an aligned `empirical` frequency attribute and a correspond
#' sample size
#'
#' @param object an `outcomes_tibble` object
#' @param i row indices or logical vector of length `nrow(object)`
#' @param j col indices, names, or logical vector of length `ncol(object)`
#' @param drop logical (default = FALSE). If `TRUE`, simplifies the result if possible.
#' @details
#' - When subsetting rows, the relative frequencies (`empirical`) are **renormalized**
#'   to sum to 1 for the selected rows.
#' - The `sampleSize` attribute is updated to reflect the sum of the original
#'   data points represented by the selected rows.
#' - Column subsetting (`j`) works like a normal tibble; attributes are preserved.
#'
#' @return An `outcomes_tibble` object with updated `empirical` and `sampleSize` attributes.
#' @noRd
#' @export
`[.outcomes_tibble` <- function(object, i, j, drop = FALSE) {
  out <- NextMethod()

  empirical <- attr(object, "empirical")[i]
  sampleSize  <- attr(object, "sampleSize")

  if (length(empirical) == 0) {
    stop("subsetting removed all rows.")
  }

  validate_outcome(
    new_outcomes_tibble(
      out,
      empirical = empirical / sum(empirical),
      sampleSize = sum(empirical * sampleSize)
    )
  )
}

#' Obtain the empirical relative frequencies
#'
#' Accessor for the `empirical` attribute of an `outcomes_tibble` object
#'
#' @param object an `outcomes_tibble` object
#' @return a numeric vector of length = `nrow(object)` with the empirical relative
#' frequencies
#'
#' @export
empirical <- function(object) UseMethod("empirical")
#' @export
empirical.outcomes_tibble <- function(object){
  attr(object, "empirical", exact = TRUE)
}

#' Set the empirical distribution
#'
#' Setter for the `empirical` attribute of an `outcomes_tibble` object
#'
#' @param object an `outcomes_tibble` object
#' @param value an non-negative numeric vector of the same length
#' as the empirical vector that sums to 1. Corresponds to the new empirical
#' distribution,
#'
#' @export
`empirical<-` <- function(object, value, tolerance = .Machine$double.eps) UseMethod("empirical<-")

#' @export
`empirical<-.outcomes_tibble` <- function(object, value, tolerance = .Machine$double.eps){
  if (!is.numeric(value)){
    stop("`value` must be a numeric vector")
  }
  if(length(attr(object, "empirical")) != length(value)){
    stop("`value` must have the same length as `empirical`")
  }
  if (any(value < 0)){
    stop("`value` must be non-negative")
  }
  if(abs(sum(value) - 1) >= sqrt(tolerance)){
    warning("`value` should sum to 1. Normalizing to 1")
    value <- value / sum(value)
  }

  attr(object, "empirical") <- value
  object
}

#' Obtain the sample size
#'
#' Accessor for the `sampleSize` attribute of an `outcomes_tibble` object
#' @param object an `outcomes_tibble` object
#' @return the sample size
#'
#' @export
sampleSize <- function(object) UseMethod("sampleSize")
#' @export
sampleSize.outcomes_tibble <- function(object){
  attr(object, "sampleSize", exact = TRUE)
}

#' Print Method for class `outcomes_tibble`
#'
#' Appends a column empirical and prints the sample size
#'
#' @param object an `outcomes_tibble` object
#' @param ... additional arguments for print
#' @return Function prints to console
#' @importFrom tibble as_tibble
#' @noRd
#' @export
print.outcomes_tibble <- function(object, ...) {
  validate_outcomes_tibble(object)
  # Create display tibble
  tbl <- as_tibble(object)
  empirical <- attr(object, "empirical", exact = TRUE)

  # Add relFreq column for display only
  tbl$empirical <- empirical

  # Use pillar for printing
  print(tbl)

  # Footer info
  cat("# Sample size:", attr(object, "sampleSize"), "\n")
  invisible(object)
}


