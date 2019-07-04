#' From MZmine output of Multiple Attributes Create Single Attribute Table
#'
#' function that takes an MZmine exported table of multiple
#' peak attributes in data.table format and separates out a table of just the
#' specified attribute with the label columns.
#'
#' @param fullAttTable data.table object; A data.table object in the format
#' exported from mzmine when multiple peak attributes are selected
#' @param firstInjColNum integer; The column index number of the first column
#' with injection peak attributes
#' @param attNum integer; the number corresponding to the position of the
#' desired attribute in the table. E.g. 3 would indicate the 3rd attribute
#' listed in the table.
#' @param numAtts integer; The number of attributes included in the table
#'
#' @import data.table
#'
#' @return Returns the single attribute table
#' @export
#'
#' @examples TODO

SingleAttTable <- function(fullAttTable,
                            firstInjColNum,
                            attNum, numAtts) {

  # column number of last label column
  lastLabelColNum <- firstInjColNum - 1

  # column number of first character column
  firstAttCol <- firstInjColNum + attNum - 1

  # total columns
  numTotAttCols <- ncol(fullAttTable)

  # attribute column sequence for specified characteristic
  attColNums <- seq(from = firstAttCol,
                     to = numTotAttCols,
                     by = numAtts)

  # create subset of imported full attribute table
  # with mz, rt, and specified attribute
  # columns for specific attribute
  dt.attTable <- subset(fullAttTable,
                         select = c(1:lastLabelColNum, attColNums))

  return(dt.attTable)

}
