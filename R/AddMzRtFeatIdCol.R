#' @title Add column to compound data table combining m/z and retention time values
#' for a unique identifier
#'
#' @description function for creating a featID column from the mz and RT
#' of each feature and placing it before the first injection
#' column
#'
#' @import data.table
#'
#' @param dt data.table object with columns corresponding the a compound's
#' mass to charge ratio (m/z) and retention time
#' @param mzColNum integer; the column index number corresponding to the m/z
#' column in \code{dt}
#' @param rtColNum integer; the column index number corresponding to the
#' retention time column in \code{dt}
#' @param placeAfterColNum integer; the new column will be inserted after the
#' column with this index. Default value is 0 putting the new column first.
#' @param colName character string; the name to be given to the new column.
#' Defaults to "mzRtFeatID"
#'
#' @return Returns the input data.table object with the column added
#' @export
#'
#' @examples TODO
AddMzRtFeatIdCol <- function(dt, mzColNum, rtColNum, placeAfterColNum = 0,
                             colName = "mzRtFeatID") {

  # make copy so original table is not altered
  dtNew <- copy(dt)

  # if the selected column name is already in the table
  if (colName %in% colnames(dtNew)) {

    stop(paste0("table already has column with name ",colName))

  }

  # final column number of new column
  finalColNumNew <- placeAfterColNum + 1

  # original number of columns in table
  ogNumCols <- ncol(dtNew)

  # convert mz and RT columns to character vectors
  mzColVect <- unlist(round(dtNew[,mzColNum, with = FALSE],4))
  rtColVect <- unlist(round(dtNew[,rtColNum, with = FALSE],3))
  mzColVect <- sprintf("%09.4f",mzColVect)
  rtColVect <- sprintf("%4.3f",rtColVect)

  # create column in data set for feature ID combining mz and RT as text
  dtNew[, (colName) :=  paste0("mz",mzColVect,"rt",rtColVect)]

  rm(mzColVect)
  rm(rtColVect)
  gc()

  # place new column in position indicated
  if (placeAfterColNum == 0) {

    setcolorder(dtNew, c(ncol(dtNew), 1:ogNumCols))

  }

  if (placeAfterColNum == ogNumCols) {
    # no order adjustment
  }

  if (placeAfterColNum != 0 & placeAfterColNum != ogNumCols) {

    colAfterNewCol <- placeAfterColNum + 1

    # Put feature ID columns before injection
    setcolorder(dtNew, c((1:placeAfterColNum),ncol(dtNew), (colAfterNewCol:ogNumCols)))
  }

  return(dtNew)

}
