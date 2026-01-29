
#' Constructor and validator for an `entities_tibble` class object
#'
#' @description
#' Creates an `entities_tibble` class object. From a `tibble` with unique rows (=entities) and
#' the empirical relative frequencies of the entities in the input `data.frame` and the sample size
#' \itemize{
#'   \item{`empirical`}{the empirical relative frequencies of the entities (unique rows)}
#'   \item{`sampleSize`}{the sample size given by the number of rows of the input `data.frame`}
#' }
#'
#' @param data A tibble-like object with unique rows
#' @param empirical a vector of empirical relative frequencies with one value per row of `data`
#' @param sampleSize a single numeric value giving the sample size. Should be at least the number of
#' rows in `data`.
#'
#' @importFrom dplyr group_by across everything summarise n select
#' @importFrom tibble new_tibble validate_tibble
#'
#' @export
new_entities_tibble <- function(data, empirical, sampleSize){
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
  validate_entities_tibble(
    new_tibble(
      data,
      nrow = nrow(data),
      empirical = empirical,
      sampleSize = sampleSize,
      class = "entities_tibble"
    )
  )
}
#' @description
#' `validate_entities_tibble()` checks an `entities_tibble` for internal consistency.
#'
#' @param object an entities_tibble class object
#' @rdname new_entities_tibble
#' @export
validate_entities_tibble <- function(object){
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
    stop("fewer observations", sampleSize, "than entities", nrow(object), call. = FALSE)
  }
  if (!is.numeric(sampleSize) || length(sampleSize) != 1) {
    stop("`sampleSize` must be a single number.")
  }

  if (sampleSize <= 0) {
    stop("`sampleSize` must be positive.")
  }

  object
}
#' Build an `entities_tibble` class object
#'
#' @description
#' `entities_tibble()` constructs a subclass of `tibble` with unique rows from an input tibble-like object. It adds two attributes
#' \itemize{
#'   \item{`empirical`}{the empirical relative frequencies of the entities (unique rows)}
#'   \item{`sampleSize`}{the sample size given by the number of rows of the input `data.frame`}
#' }
#'
#' @param data a tibble-like object
#'
#'
#'
#' @export
entities_tibble <- function(data){
  if (!is.data.frame(data)){
    stop("df should be a data frame but is", class(data))
  }
  ## the number of rows is the sample size
  sampleSize <- nrow(data)
  ## group data by all columns
  data <- group_by(data, across(everything()))
  ## summarize and drop
  data <- summarise(data, empirical = n() / sampleSize, .groups = "drop")

  new_entities_tibble(
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
#' @param object an `entities_tibble` object
#' @param i row indices selected by an upstream `dplyr` function
#' @param ... additional parameters
#'
#' @importFrom dplyr dplyr_row_slice
#' @noRd
#' @export

dplyr_row_slice.entities_tibble <- function(object, i, ...) {
  if (length(i) == 0) {
    stop("subsetting removed all rows.")
  }

  out <- NextMethod()

  empirical <- attr(object, "empirical")[i]
  sampleSize  <- attr(object, "sampleSize")


  validate_entities(
    new_entities(
      out,
      empirical = empirical / sum(empirical),
      sampleSize = sum(empirical * sampleSize)
    )
  )
}

#' Subset an `entities_tibble`
#'
#' Subsetting an `entities_tibble` works similarily to subsetting a `tibble`,
#' but the object maintains an aligned `empirical` frequency attribute and a correspond
#' sample size
#'
#' @param object an `entities_tibble` object
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
#' @return An `entities_tibble` object with updated `empirical` and `sampleSize` attributes.
#'#' @export
`[.entities_tibble` <- function(object, i, j, drop = FALSE) {
  out <- NextMethod()

  empirical <- attr(object, "empirical")[i]
  sampleSize  <- attr(object, "sampleSize")

  if (length(empirical) == 0) {
    stop("subsetting removed all rows.")
  }

  validate_entities(
    new_entities(
      out,
      empirical = empirical / sum(empirical),
      sampleSize = sum(empirical * sampleSize)
    )
  )
}

#' Obtain the empirical relative frequencies
#'
#' Accessor for the `empirical` attribute of an `entities_tibble` object
#'
#' @param object an `entities_tibble` object
#' @return a numeric vector of length = `nrow(object)` with the empirical relative
#' frequencies
#'
#' @export
empirical <- function(object) UseMethod("empirical")
#' @export
empirical.entities_tibble <- function(object){
  attr(object, "empirical", exact = TRUE)
}

#' Obtain the sample size
#'
#' Accessor for the `sampleSize` attribute of an `entities_tibble` object
#' @param object an `entities_tibble` object
#' @return the sample size
#'
#' @export
sampleSize <- function(object) UseMethod("sampleSize")
#' @export
sampleSize.entities_tibble <- function(object){
  attr(object, "sampleSize", exact = TRUE)
}

#' Print Method for class `entities_tibble`
#'
#' Appends a column empirical and prints the sample size
#'
#' @param object an `entities_tibble` object
#' @param ... additional arguments for print
#' @return Function prints to console
#' @export
print.entities_tibble <- function(object, ...) {
  validate_entities(object)
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


