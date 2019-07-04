#' Create PDF file from Plot List
#'
#' This function creates a pdf of the plots in list.plotPage
#' and exports them to the file path/name PdfFilePath
#'
#' @param list.plotPage list of lists of ggplot plot objects
#' @param orientation "portrait" or "landscape" to indicate orientation of outuput
#' pdf
#' @param folderPath path to directory where new pdf file will be generated
#' @param name name of output PDF file. ".pdf" suffix will be added
#'
#' @return none
#' @export
#'
#' @examples TODO
PdfPlotsSingleFile <- function(list.plotPage, orientation, folderPath, name) {

  # list.plotPage is the list of lists of plot objects
  # orientation is the orientation of the paper
  # "portrait" or "landscape"
  # folderPath = path of folder to write to
  # name = base name of output file, time stamp
  # and .csv are added

  PdfFilePath <- paste0(folderPath, "/", name, ".pdf")

  # paper orientation
  if (orientation == "portrait") {
    paper <- "letter"
    width <- 7.5
    height <- 10
  } else if (orientation == "landscape") {
    paper <- "USr"
    width <- 10
    height <- 7.5
  } else {
    stop("please enter portrait or landscape orientation")
  }

  pdf(PdfFilePath, width = width, height = height,
      paper = paper, onefile = TRUE)
  for (i in seq(length(list.plotPage))) {
    gridExtra::grid.arrange(grobs = list.plotPage[[i]], ncol = 1)

    percComplete <- i/length(list.plotPage) * 100

    print(paste0(name, " pdf completion ",
                 format(percComplete, nsmall = 1),"%"))
  }
  dev.off()
}
