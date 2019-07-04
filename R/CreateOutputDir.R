#' Create a new directory for output files
#'
#' @param dirName the name of the directory to be created
#' @param parentDir the directory (folder) path where the output folder
#' is to be created
#'
#' @return returns the path of the newly created directory
#' @export


CreateOutputDir <- function(parentDir, dirName) {

  startTimeStamp <- format(Sys.time(), "(%Y.%m.%d.%H%M)")

  OutputFolderName <- paste0(dirName,startTimeStamp)

  # create output folder if it doesn't exist
  outputFolderPath <- paste0(parentDir,"/",OutputFolderName)

  if (!dir.exists(outputFolderPath)) {
    dir.create(outputFolderPath)
  }

  return(outputFolderPath)

}

