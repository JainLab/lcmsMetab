#' Get runtime of script
#'
#' Get the system time that has elapsed
#' since the system time of the startTime
#' argument
#' @param startTime sys.time() generated value from when timing should start
#' @return The elapsed system time since startTime
#' @export

runTime <- function(startTime){
  format(Sys.time() -  startTime, format = "")
}
