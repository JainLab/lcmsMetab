#' Write Data Table to CSV file
#'
#' @param table data table object
#' @param folderPath path to directory where table will be written
#' @param name name of file to be written - ".csv" will be appended
#'
#' @return Does not retun value. Writes csv file.
#' @export
#'
#' @examples TODO
WrtTable <- function(table, folderPath, name) {

  # export table
  # add data and file extension to the file name
  OutputFileName <- paste0(folderPath,"/",
                           name,".csv")

  #Write the adjusted dataframe to csv
  data.table::fwrite(table,
         file = OutputFileName,
         row.names = FALSE)

}
