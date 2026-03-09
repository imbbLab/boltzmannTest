#' Construct for an `outcomes_tibble` object
#'
#' Constructs an `outcomes_tibble`  object with unique
#' rows from an input tibble-like object. It adds two attributes
#' \itemize{
#'   \item{`empirical` the empirical relative frequencies of the outcomes
#'   (unique rows)}
#'   \item{`sampleSize` the sample size given by the number of rows of the
#'   input `data.frame`}
#' }
#'
#' @param data a data frame with outcomes in the rows
#' @param empirical a numeric vector with one value per row of `data`
#' @param sampleSize an integer value indicating the sample Size
#'
#' @returns An `outcomes_tibble` object
#'
#' @importFrom dplyr group_by across everything summarise n select
#' @importFrom tibble new_tibble validate_tibble
#'
#' @keywords internal
new_outcomes_tibble <- function(data, empirical, sampleSize){

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
#' Validate an `outcomes_tibble` class object
#'
#' Checks an `outcomes_tibble` for internal consistency.
#'
#' @param object an outcomes_tibble class object
#'
#' @returns An `outcomes_tibble` object. If validation fails it stops with
#' an error message about what could not be validated.
#'
#' @keywords internal
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
    stop("`empirical` must have one value per row of `data`. There are ",
         length(empirical), " entries in `empirical` and ", nrow(object),
         " rows in `data`")
  }
  if (any(empirical < 0)){
    stop("`empirical` contains negative entries")
  }
  if (abs(sum(empirical) - 1) > sqrt(.Machine$double.eps)){
    stop("`empirical` does not sum to 1")
  }
  if (!is.numeric(sampleSize) || length(sampleSize) != 1) {
    stop("`sampleSize` must be a single number.")
  }
  if (sampleSize <= 0) {
    stop("`sampleSize` must be positive.")
  }
  if (sampleSize < nrow(object)){
    stop("fewer observations (", sampleSize,
         ") than outcomes (", nrow(object), ")")
  }


  object
}

#' Build an `outcomes_tibble` object
#'
#' `outcomes_tibble()` constructs an `outcomes_tibble()` from a data.frame-like
#' object. It will have two attributes:
#'  * `empirical` the empirical relative frequencies of the outcomes
#'    (unique rows)
#'  * `sampleSize` the sample size given by the number of rows of the
#'    input `data.frame`
#'
#' @param data a data frame
#'
#' @examples
#' ## use the data from NHANES
#' data(nhanes)
#' ## select only BMXWT and RIAGENDR columns
#' data <- subset(nhanes, select = c(BMXWT, RIAGENDR))
#' ## look at data
#' str(data)
#' ## create an outcomes_tibble object
#' outcomes <- outcomes_tibble(data)
#' str(outcomes)
#'
#'
#'
#' @returns An `outcomes_tibble` object
#'
#' @export
outcomes_tibble <- function(data){
  ## check if data is a data.frame
  if (!is.data.frame(data)){
    stop("df should be a data frame but is", class(data))
  }
  ## detect and remove cases with non-finite and NA entries
  completeCases <- completeCases_internal(data)
  if (any(! completeCases)){
    warning("there are rows with non-finite and/or NA values. Removing them.")
    data <- data[completeCases, ]

  }
  ## convert character variables to factors
  data <- data.frame(
    as.list(data),
    stringsAsFactors = TRUE,
    check.names = FALSE
  )
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

#' `dplyr` row slicing for an `outcomes_tibble`
#'
#'
#' This function is automatically called after the `dplyr` functions `filter`,
#' `slice`, `arrange` etc. It automatically updates (marginalize) the empirical
#' relative frequencies and the sample size
#'
#' @param object an `outcomes_tibble` object
#' @param i row indices selected by an upstream `dplyr` function
#' @param ... additional parameters
#'
#' @returns A sliced `outcomes_tibble`
#'
#' @importFrom dplyr dplyr_row_slice
#' @export
dplyr_row_slice.outcomes_tibble <- function(object, i, ...) {

  if (length(i) == 0 || all(i == FALSE)) {
    stop("subsetting removed all rows.")
  }

  out <- NextMethod()

  empirical <- attr(object, "empirical")[i]
  sampleSize  <- attr(object, "sampleSize")


  validate_outcomes_tibble(
    new_outcomes_tibble(
      out,
      empirical = empirical / sum(empirical),
      sampleSize = sum(empirical * sampleSize)
    )
  )
}

#' Subset an `outcomes_tibble`
#'
#' Subsetting an `outcomes_tibble` works similarily to subsetting a `tibble`,
#' but the object maintains an aligned `empirical` frequency attribute and
#' sample size
#'
#' @param object an `outcomes_tibble` object
#' @param i row indices or logical vector of length `nrow(object)`
#' @param j col indices, names, or logical vector of length `ncol(object)`
#' @param drop logical (default = FALSE). If `TRUE`, simplifies the result if
#' possible.
#' @details
#' - When subsetting rows, the relative frequencies (`empirical`) are
#'   **renormalized**
#'   to sum to 1 for the selected rows.
#' - The `sampleSize` attribute is updated to reflect the sum of the original
#'   data points represented by the selected rows.
#' - Column subsetting (`j`) works like a normal tibble; attributes are
#'   preserved.
#'
#' @return An `outcomes_tibble` object with updated `empirical` and
#' `sampleSize` attributes.
#'
#' @export
`[.outcomes_tibble` <- function(object, i, j, drop = FALSE) {
  out <- NextMethod()

  empirical <- attr(object, "empirical")[i]
  sampleSize  <- attr(object, "sampleSize")

  if (length(empirical) == 0) {
    stop("subsetting removed all rows.")
  }


  ## group data by all remaining columns
  out <- group_by(out, across(everything()))
  empirical <- vapply(
    X = attr(out, "groups")$.rows,
    FUN = function(i){
      sum(empirical[i])
    },
    FUN.VALUE = 1
  )

  out <- attr(out, "groups")
  out$.rows <- NULL
  attr(out, ".drop") <- NULL
  validate_outcomes_tibble(
    new_outcomes_tibble(
      out,
      empirical = empirical / sum(empirical),
      sampleSize = sum(empirical * sampleSize)
    )
  )
}

#' Get or set the `empirical` attribute of an `outcomes_tibble` object
#'
#' @description
#' `empirical()` returns the `empirical` attribute of an `outcomes_tibble`
#' object, i.e. the empirical frequencies of the outcomes
#'
#' @param object an `outcomes_tibble` object
#' @return a numeric vector of length = `nrow(object)` with the empirical
#' relative frequencies
#'
#' @rdname empirical
#' @export
empirical <- function(object) UseMethod("empirical")
#' @export
empirical.outcomes_tibble <- function(object){
  attr(object, "empirical", exact = TRUE)
}

#' empirical<-
#'
#' @description
#' `empirical<-` sets the `empirical` attribute of an `outcomes_tibble` object
#'
#' @param object an `outcomes_tibble` object
#' @param value an non-negative numeric vector of the same length
#' as the empirical vector that sums to 1. Corresponds to the new empirical
#' distribution,
#'
#' @rdname empirical
#' @export
`empirical<-` <- function(object, value) UseMethod("empirical<-")

#' @export
`empirical<-.outcomes_tibble` <- function(object, value){
  if (!is.numeric(value)){
    stop("`value` must be a numeric vector")
  }
  if(length(attr(object, "empirical")) != length(value)){
    stop("`value` must have the same length as `empirical`")
  }
  if (any(value < 0)){
    stop("`value` must be non-negative")
  }
  if(abs(sum(value) - 1) >= sqrt(.Machine$double.eps)){
    warning("`value` should sum to 1. Normalizing to 1")
    value <- value / sum(value)
  }

  attr(object, "empirical") <- value
  object
}

#' Get the sample size
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
#'
#' @examples
#'
#' data(kidneyStones)
#' outcomes <- outcomes_tibble(kidneyStones)
#' print(outcomes)
#'
#'
#' @importFrom tibble as_tibble
#'
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

completeCases_internal <- function(data){
  ## check for rows with NAs
  completeCases <- do.call(
    "rbind",
    lapply(
      data,
      function(x){
        is.finite(x) | (! is.numeric(x) & !is.na(x))
      }
    )
  )
  apply(completeCases, 2, all)
}

