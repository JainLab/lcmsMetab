#' Replace 0 Cells and Empty Cells in data.table with NA Remove Duplicate Rows
#'
#' @param x data table object
#'
#' @import data.table
#'
#' @return returns the revised data table object
#' @export
#'
#' @examples TODO

NaZeroEmptyRmDup <- function(x) {

  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package \"data.table\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if (data.table::is.data.table(x) == FALSE) {
    stop("input must be data.table object")

  }

  # input x is data table

  # change 0 cells to NA
  x[x == 0] <- NA
  # change empty cells to NA
  x[x == ""] <- NA

  # remove any columns that are all NA - resulting from blank column
  x <- x[,which(unlist(lapply(x, function(x)!all(is.na(x))))), with = FALSE]
  # x <- x[,which(unlist(lapply(x, function(x)!all(is.na(x)))))]

  # remove duplicate rows
  x <- unique(x)

  return(x)

}
