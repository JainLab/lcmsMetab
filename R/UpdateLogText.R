#' Update summary text file
#'
#' Add text to log file and update text log file - used mostly for development
#'
#' @param newText new text to add to the log file
#' @param ... additional lines of text
#' @param logText character vector of existing text of log file
#'
#' @return Returns text for the log file
#' @export
#'
#' @examples summaryText <-
#' updateSummary("Initialize Summary File")


UpdateLogText <- function(logText, newText, ...){

  # create timestamp for when new text is added
  timeStamp <- format(Sys.time(), "%Y.%m.%d.%H%M")

  # add new text to existing logText and add
  # timeStamp
  logText <- c(logText,"",newText,..., timeStamp)

  return(logText)
}
