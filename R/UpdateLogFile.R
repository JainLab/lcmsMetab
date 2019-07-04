#' create log text file and Update with new text
#'
#' @param logFilePath file path of log file
#' @param logText text of log file
#'
#' @return No return - simply updates log file
#' @export
#'
#' @examples
#'
UpdateLogFile <- function(logFilePath, logText){
  logFile <- file(logFilePath)
  writeLines(logText, logFile)
  close(logFile)
}
