#' Converts data.table Columns to Z-scores by Column
#'
#' converts values to z scores by column for all columns between the column
#' indices indicated
#'
#' @param dt data.table object
#' @param indFirstConvCol numeric; column index number of first column to be converted
#' z-scores. Defaults to 1.
#' @param indLastConvCol numeric; column index number of last column to be converted.
#' Defaults to the last column in the data table.
#'
#' @import data.table
#'
#' @return Returns the adjusted data.table object
#' @export
#'
#' @examples TODO

ZScoreCols <- function(dt, indFirstConvCol = 1, indLastConvCol = ncol(dt)) {

  for (colj in indFirstConvCol:indLastConvCol) {

    set(dt,
        i = NULL, # all rows
        j = colj, # iterated column
        value = scale(dt[[colj]], center = TRUE, scale = TRUE))
  }

  return(dt)

}
