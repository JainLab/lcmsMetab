#' Create Peak Prediction Table
#'
#' This function creates a table of peak group statistics used to generate
#' models of peak quality and to predict peak quality from those models.
#'
#' @param dt.AllAtt data.table object; contains columns:  row m/z
#' row retention time, for the data file elements select:
#' peak duration,
#' peak height,
#' peak area,
#' peak FWHM,
#' peak tailing factor,
#' peak asymmetry factor
#' @param mzColNum column index number of column containing the compound mass
#' to charge ratio. Defaults to 1.
#' @param rtColNum column index number of column containing retention time of
#' the compound. Defaults to 2.
#'
#' @return returns a data.table object of prediction parameters
#' @export
#'
#' @examples TODO

PeakPredictionTable <- function(dt.AllAtt, mzColNum = 1, rtColNum = 2) {

  # function timing --------------------------------------

  startTime <- Sys.time()

  # replace zeros and empty cells, remove empty columns -------------

  # apply defined function
  dt.AllAtt <- NaZeroEmptyRmDup(dt.AllAtt)

  # if columns were imported as non-numeric convert to numeric --------

  # in certain cases it seems the wrong class is applied to columns
  # during import because of the way that fread functions, sparse
  # numeric columns are not identified as numeric

  # convert all columns to numeric
  for (col in colnames(dt.AllAtt)) {
    set(dt.AllAtt, j = col, value = as.numeric(dt.AllAtt[[col]]))
  }

  print("replace 0 and empty values, correct column classes")
  print(runTime(startTime))

  # create individual characteristic tables -----------------------------

  # duration, height, area, FWHM, tailing factor,
  # asymmetry factor
  numChars <- 6

  # use defined function to create individual
  # characteristic tables
  dt.PkDur <- SingleAttTable(dt.AllAtt, firstInjColNum,
                             1, numChars)
  dt.PkHt <- SingleAttTable(dt.AllAtt, firstInjColNum,
                            2, numChars)
  dt.PkArea <- SingleAttTable(dt.AllAtt, firstInjColNum,
                              3, numChars)
  dt.PkFWHM <- SingleAttTable(dt.AllAtt, firstInjColNum,
                              4, numChars)
  dt.PkTfact <- SingleAttTable(dt.AllAtt, firstInjColNum,
                               5, numChars)
  dt.PkAfact <- SingleAttTable(dt.AllAtt, firstInjColNum,
                               6, numChars)

  print("separate attributes")
  print(runTime(startTime))

  # add featID column to each table -----------------------------------------



  # apply addFeatID function defined above to tables
  dt.PkHt <- AddMzRtFeatIdCol(dt.PkHt, mzColNum, rtColNum,
                              placeAfterColNum = 0, colName = "featID")
  dt.PkArea <- AddMzRtFeatIdCol(dt.PkArea, mzColNum, rtColNum,
                                placeAfterColNum = 0, colName = "featID")
  dt.PkDur <- AddMzRtFeatIdCol(dt.PkDur, mzColNum, rtColNum,
                               placeAfterColNum = 0, colName = "featID")
  dt.PkFWHM <- AddMzRtFeatIdCol(dt.PkFWHM, mzColNum, rtColNum,
                                placeAfterColNum = 0, colName = "featID")
  dt.PkTfact <- AddMzRtFeatIdCol(dt.PkTfact, mzColNum, rtColNum,
                                 placeAfterColNum = 0, colName = "featID")
  dt.PkAfact <- AddMzRtFeatIdCol(dt.PkAfact, mzColNum, rtColNum,
                                 placeAfterColNum = 0, colName = "featID")

  # create feature summary table ----------------------------------

  # update column numbers
  mzColNumNew <- mzColNum + 1

  rtColNumNew <- rtColNum + 1

  # take relevant column from dt.PkHt
  dt.featSummary <- subset(dt.PkHt,
                           select = c(1, mzColNumNew, rtColNumNew))

  # convert long format --------------------------------------------

  # apply to tables
  # convert to long form with melt
  firstInjNewColNum <- firstInjColNum + 1

  dt.PkHt.long <- melt(dt.PkHt, id.vars = "featID",
                       measure.vars = c(firstInjNewColNum:ncol(dt.PkHt)),
                       variable.name = "injection",
                       value.name = "pkHt",
                       na.rm = TRUE)

  dt.PkArea.long <- melt(dt.PkArea, id.vars = "featID",
                         measure.vars = c(firstInjNewColNum:ncol(dt.PkArea)),
                         variable.name = "injection",
                         value.name = "pkArea",
                         na.rm = TRUE)

  dt.PkDur.long <- melt(dt.PkDur, id.vars = "featID",
                        measure.vars = c(firstInjNewColNum:ncol(dt.PkDur)),
                        variable.name = "injection",
                        value.name = "pkDur",
                        na.rm = TRUE)

  dt.PkFWHM.long <- melt(dt.PkFWHM, id.vars = "featID",
                         measure.vars = c(firstInjNewColNum:ncol(dt.PkFWHM)),
                         variable.name = "injection",
                         value.name = "pkFWHM",
                         na.rm = TRUE)

  dt.PkTfact.long <- melt(dt.PkTfact, id.vars = "featID",
                          measure.vars = c(firstInjNewColNum:ncol(dt.PkTfact)),
                          variable.name = "injection",
                          value.name = "pkTfact",
                          na.rm = TRUE)

  dt.PkAfact.long <- melt(dt.PkAfact, id.vars = "featID",
                          measure.vars = c(firstInjNewColNum:ncol(dt.PkAfact)),
                          variable.name = "injection",
                          value.name = "pkAfact",
                          na.rm = TRUE)

  print("long format tables generated")
  print(runTime(startTime))

  # remove injection column suffix ------------------------------

  dt.PkHt.long[, injection := gsub("mzXML.*$","mzXML",injection)]
  dt.PkArea.long[, injection := gsub("mzXML.*$","mzXML",injection)]
  dt.PkDur.long[, injection := gsub("mzXML.*$","mzXML",injection)]
  dt.PkFWHM.long[, injection := gsub("mzXML.*$","mzXML",injection)]
  dt.PkTfact.long[, injection := gsub("mzXML.*$","mzXML",injection)]
  dt.PkAfact.long[, injection := gsub("mzXML.*$","mzXML",injection)]

  # find total num injections (with features) ------------

  numTotInj <- uniqueN(dt.PkHt.long, by = c("injection"))

  # proportion detected ----------------------------------

  # create table of hits per featID
  dt.featHits <- dt.PkHt.long[,sum(!is.na(pkHt)),
                              keyby = .(featID)]
  colnames(dt.featHits)[2] <- "numHits"

  # create column for proportion detected
  dt.featHits[, propDetected := numHits / numTotInj]

  # add to dt.featSummary
  setkey(dt.featSummary, featID)
  dt.featSummary[dt.featHits, propDetected := propDetected]

  rm(dt.featHits)
  gc()

  # height -------------------------------

  # create table for max height
  dt.maxHeight <- dt.PkHt.long[,max(pkHt, na.rm = TRUE),
                               keyby = .(featID)]
  colnames(dt.maxHeight)[2] <- "maxHeight"

  # get log maxHeight
  dt.maxHeight[, logMaxHt := log(maxHeight)]

  # create table for median height
  dt.medHeight <- dt.PkHt.long[,median(pkHt, na.rm = TRUE),
                               keyby = .(featID)]
  colnames(dt.medHeight)[2] <- "medHeight"

  # get log medHeight
  dt.medHeight[, logMedHt := log(medHeight)]

  # table for standard deviation of height
  dt.sdHeight <- dt.PkHt.long[, sd(pkHt, na.rm = TRUE),
                              keyby = .(featID)]
  colnames(dt.sdHeight)[2] <- "sdHeight"

  # log sd height
  dt.sdHeight[, logSdHeight := log(sdHeight)]

  # add to dt.featSummary
  setkey(dt.featSummary, featID)
  dt.featSummary[dt.maxHeight, logMaxHt := logMaxHt]
  dt.featSummary[dt.medHeight, logMedHt := logMedHt]
  dt.featSummary[dt.sdHeight, logSdHeight := logSdHeight]

  rm(dt.maxHeight)
  rm(dt.medHeight)
  gc()

  # area -------------------------------

  # create table for max area
  dt.maxArea <- dt.PkArea.long[,max(pkArea, na.rm = TRUE),
                               keyby = .(featID)]
  colnames(dt.maxArea)[2] <- "maxArea"

  # get log maxArea
  dt.maxArea[, logMaxArea := log(maxArea)]

  # create table for median area
  dt.medArea <- dt.PkArea.long[,median(pkArea, na.rm = TRUE),
                               keyby = .(featID)]
  colnames(dt.medArea)[2] <- "medArea"

  # get log medArea
  dt.medArea[, logMedArea := log(medArea)]

  # table for standard deviation of height
  dt.sdArea <- dt.PkArea.long[, sd(pkArea, na.rm = TRUE),
                              keyby = .(featID)]
  colnames(dt.sdArea)[2] <- "sdArea"

  # log sd height
  dt.sdArea[, logSdArea := log(sdArea)]

  # add to dt.featSummary
  setkey(dt.featSummary, featID)
  dt.featSummary[dt.maxArea, logMaxArea := logMaxArea]
  dt.featSummary[dt.medArea, logMedArea := logMedArea]
  dt.featSummary[dt.sdArea, logSdArea := logSdArea]

  rm(dt.maxArea)
  rm(dt.medArea)
  gc()

  # duration ----------------------------

  # replace durations longer than max allowable with max allowable
  dt.PkDur.long[, pkDur := ifelse(pkDur > maxDurSet, maxDurSet, pkDur)]

  # create table for max duration
  dt.maxDur <- dt.PkDur.long[,max(pkDur, na.rm = TRUE),
                             keyby = .(featID)]
  colnames(dt.maxDur)[2] <- "maxDur"

  # get log maxDur
  dt.maxDur[, logMaxDur := log(maxDur)]

  # create table for median duration
  dt.medDur <- dt.PkDur.long[,median(pkDur, na.rm = TRUE),
                             keyby = .(featID)]
  colnames(dt.medDur)[2] <- "medDur"

  # get log maxDur
  dt.medDur[, logMedDur := log(medDur)]

  # create table for mean duration
  dt.avgDur <- dt.PkDur.long[,mean(pkDur, na.rm = TRUE),
                             keyby = .(featID)]
  colnames(dt.avgDur)[2] <- "avgDur"

  # create table for standard deviation of duration
  dt.sdDur <- dt.PkDur.long[,sd(pkDur, na.rm = TRUE),
                            keyby = .(featID)]
  colnames(dt.sdDur)[2] <- "stDevDur"

  # to prevent issues with other calculations replace stDevDur = 0
  # with very low value
  dt.sdDur[,stDevDur := ifelse(stDevDur == 0, 0.000001, stDevDur)]

  # get log stDevDur
  dt.sdDur[, logStDevDur := log(stDevDur)]

  # add avgDur to dt.sdDur
  dt.sdDur[dt.avgDur, avgDur := avgDur]

  # create column for coefficient of variation
  dt.sdDur[, CvDur := stDevDur / avgDur]

  # add maxDur to dt.sdDur
  dt.sdDur[dt.maxDur, maxDur := maxDur]

  # create column for stDevDur / maxDur
  dt.sdDur[, custVarMeasDur := stDevDur / maxDur]

  # add to dt.featSummary
  setkey(dt.featSummary, featID)
  dt.featSummary[dt.maxDur, logMaxDur := logMaxDur]
  dt.featSummary[dt.medDur, logMedDur := logMedDur]
  dt.featSummary[dt.sdDur, logStDevDur := logStDevDur]
  dt.featSummary[dt.sdDur, CvDur := CvDur]
  dt.featSummary[dt.sdDur, custVarMeasDur := custVarMeasDur]

  rm(dt.maxDur)
  rm(dt.medDur)
  gc()

  # FWHM ----------------------------

  # replace FWHM longer than max allowable with max allowable
  dt.PkFWHM.long[, pkFWHM := ifelse(pkFWHM > maxDurSet, maxDurSet, pkFWHM)]

  # create table for max FWHM
  dt.maxFWHM <- dt.PkFWHM.long[,max(pkFWHM, na.rm = TRUE),
                               keyby = .(featID)]
  colnames(dt.maxFWHM)[2] <- "maxFWHM"

  # get log maxFWHM
  dt.maxFWHM[, logMaxFWHM := log(maxFWHM)]

  # create table for median FWHM
  dt.medFWHM <- dt.PkFWHM.long[,median(pkFWHM, na.rm = TRUE),
                               keyby = .(featID)]
  colnames(dt.medFWHM)[2] <- "medFWHM"

  # get log medFWHM
  dt.medFWHM[, logMedFWHM := log(medFWHM)]

  # add to dt.featSummary
  setkey(dt.featSummary, featID)
  dt.featSummary[dt.maxFWHM, logMaxFWHM := logMaxFWHM]
  dt.featSummary[dt.medFWHM, logMedFWHM := logMedFWHM]

  rm(dt.maxFWHM)
  rm(dt.medFWHM)
  gc()

  # area/height comparisons -------------------------------

  # add peak height to dt.PkArea.long
  setkeyv(dt.PkHt.long, c("featID", "injection"))
  setkeyv(dt.PkArea.long, c("featID", "injection"))
  dt.PkArea.long[dt.PkHt.long, pkHt := pkHt]

  # create column for area/height
  dt.PkArea.long[, relArea := pkArea / pkHt]

  # create table for max relArea
  dt.maxRelArea <- dt.PkArea.long[,max(relArea, na.rm = TRUE),
                                  keyby = .(featID)]
  colnames(dt.maxRelArea)[2] <- "maxRelArea"

  # get log maxRelArea
  dt.maxRelArea[, logMaxRelArea := log(maxRelArea)]

  # create table for median relArea
  dt.medRelArea <- dt.PkArea.long[,median(relArea, na.rm = TRUE),
                                  keyby = .(featID)]
  colnames(dt.medRelArea)[2] <- "medRelArea"

  # get log medRelArea
  dt.medRelArea[, logMedRelArea := log(medRelArea)]

  # table for min relArea
  dt.minRelArea <- dt.PkArea.long[,log(min(relArea, na.rm = TRUE)),
                                  keyby = .(featID)]
  colnames(dt.minRelArea)[2] <- "logMinRelArea"

  # table for mean relArea
  dt.avgRelArea <- dt.PkArea.long[,mean(relArea, na.rm = TRUE),
                                  keyby = .(featID)]
  colnames(dt.avgRelArea)[2] <- "avgRelArea"

  # table for standard deviation of relArea
  dt.sdRelArea <- dt.PkArea.long[,sd(relArea, na.rm = TRUE),
                                 keyby = .(featID)]
  colnames(dt.sdRelArea)[2] <- "sdRelArea"

  # get log of sdRelArea
  dt.sdRelArea[, logSdRelArea := log(sdRelArea)]

  # add average column to standard deviation table
  dt.sdRelArea[dt.avgRelArea, avgRelArea := avgRelArea]

  # create coefficient of variation
  dt.sdRelArea[, logCvRelArea := log(sdRelArea/avgRelArea)]

  # add maxRelArea to dt.sdRelArea
  dt.sdRelArea[dt.maxRelArea, maxRelArea := maxRelArea]

  # create column for sdRelArea / maxRelArea
  dt.sdRelArea[, custVarMeasRelArea := sdRelArea / maxRelArea]

  # add to dt.featSummary
  setkey(dt.featSummary, featID)
  dt.featSummary[dt.maxRelArea, logMaxRelArea := logMaxRelArea]
  dt.featSummary[dt.medRelArea, logMedRelArea := logMedRelArea]
  dt.featSummary[dt.minRelArea, logMinRelArea := logMinRelArea]
  dt.featSummary[dt.sdRelArea, logSdRelArea := logSdRelArea]
  dt.featSummary[dt.sdRelArea, logCvRelArea := logCvRelArea]
  dt.featSummary[dt.sdRelArea, custVarMeasRelArea := custVarMeasRelArea]

  rm(dt.maxRelArea)
  rm(dt.medRelArea)
  gc()

  # sdArea / sdHeight ----------------------------------------

  # add sdArea to dt.sdHeight
  dt.sdHeight[dt.sdArea, sdArea := sdArea]

  # create column for sdArea/sdHeight
  dt.sdHeight[, sdRatioAreaHt := sdHeight / sdArea]

  # add to dt.featSummary
  setkey(dt.featSummary, featID)
  dt.featSummary[dt.sdHeight, sdRatioAreaHt := sdRatioAreaHt]

  rm(dt.sdHeight)
  rm(dt.sdArea)
  gc()

  # FWHM/Duration comparisons -------------------------------

  # add FWHM to dt.PkDur.Long
  setkeyv(dt.PkDur.long, c("featID", "injection"))
  setkeyv(dt.PkFWHM.long, c("featID", "injection"))
  dt.PkDur.long[dt.PkFWHM.long, pkFWHM := pkFWHM]

  # create column for area/height
  dt.PkDur.long[, relFWHM := pkFWHM / pkDur]

  # create table for max relFWHM
  dt.maxRelFWHM <- dt.PkDur.long[,max(relFWHM, na.rm = TRUE),
                                 keyby = .(featID)]
  colnames(dt.maxRelFWHM)[2] <- "maxRelFWHM"

  # get log maxRelFWHM
  dt.maxRelFWHM[, logMaxRelFWHM := log(maxRelFWHM)]

  # create table for median relFWHM
  dt.medRelFWHM <- dt.PkDur.long[,median(relFWHM, na.rm = TRUE),
                                 keyby = .(featID)]
  colnames(dt.medRelFWHM)[2] <- "medRelFWHM"

  # add column for log(medRelFHWM)
  dt.medRelFWHM[, logMedRelFWHM := log(medRelFWHM)]

  # create table for minimum relFWHM
  dt.minRelFWHM <- dt.PkDur.long[,min(relFWHM, na.rm = TRUE),
                                 keyby = .(featID)]
  colnames(dt.minRelFWHM)[2] <- "minRelFWHM"

  # create table for mean relFWHM
  dt.avgRelFWHM <- dt.PkDur.long[,mean(relFWHM, na.rm = TRUE),
                                 keyby = .(featID)]
  colnames(dt.avgRelFWHM)[2] <- "avgRelFWHM"

  # create table for standard deviation of relFWHM
  dt.sdRelFWHM <- dt.PkDur.long[,sd(relFWHM, na.rm = TRUE),
                                keyby = .(featID)]
  colnames(dt.sdRelFWHM)[2] <- "sdRelFWHM"

  # to prevent issues with other calculations replace sdRelFWHM = 0
  # with very low value
  dt.sdRelFWHM[,sdRelFWHM := ifelse(sdRelFWHM == 0, 0.000001, sdRelFWHM)]

  # get log sd RelFWHM
  dt.sdRelFWHM[,logSdRelFWHM := log(sdRelFWHM)]

  # add mean to standard deviation table
  dt.sdRelFWHM[dt.avgRelFWHM, avgRelFWHM := avgRelFWHM]

  # create coefficient of variation column
  dt.sdRelFWHM[,CvRelFWHM := sdRelFWHM/avgRelFWHM]

  # create log cv column
  dt.sdRelFWHM[,logCvRelFWHM := log(CvRelFWHM)]

  # add logMaxRelFWHM to dt.sdRelFWHM
  dt.sdRelFWHM[dt.maxRelFWHM, maxRelFWHM := maxRelFWHM]

  # create custom relative variance measure
  dt.sdRelFWHM[,custVarMeasRelFWHM := sdRelFWHM / maxRelFWHM]

  # add to dt.featSummary
  setkey(dt.featSummary, featID)
  dt.featSummary[dt.maxRelFWHM, logMaxRelFWHM := logMaxRelFWHM]
  dt.featSummary[dt.medRelFWHM, medRelFWHM := medRelFWHM]
  dt.featSummary[dt.medRelFWHM, logMedRelFWHM := logMedRelFWHM]
  dt.featSummary[dt.minRelFWHM, minRelFWHM := minRelFWHM]
  dt.featSummary[dt.sdRelFWHM, logSdRelFWHM := logSdRelFWHM]
  dt.featSummary[dt.sdRelFWHM, logCvRelFWHM := logCvRelFWHM]
  dt.featSummary[dt.sdRelFWHM, custVarMeasRelFWHM := custVarMeasRelFWHM]

  rm(dt.maxRelFWHM)
  rm(dt.medRelFWHM)
  rm(dt.avgRelFWHM)
  rm(dt.sdRelFWHM)
  gc()

  # area/FWHM (specWidth) ---------------------------------------------

  # add FWHM to dt.PkArea.Long
  setkeyv(dt.PkArea.long, c("featID", "injection"))
  setkeyv(dt.PkFWHM.long, c("featID", "injection"))
  dt.PkArea.long[dt.PkFWHM.long, pkFWHM := pkFWHM]

  # create area/FWHM
  dt.PkArea.long[, specWidth := pkFWHM / pkArea]

  # create table for standard deviation of relFWHM
  dt.sdSpecWidth <- dt.PkArea.long[,sd(specWidth, na.rm = TRUE),
                                   keyby = .(featID)]
  colnames(dt.sdSpecWidth)[2] <- "sdSpecWidth"

  # create column for log of sdSpecWidth
  dt.sdSpecWidth[,logSdSpecWidth := log(sdSpecWidth)]

  # add to dt.featSummary
  setkey(dt.featSummary, featID)
  dt.featSummary[dt.sdSpecWidth, logSdSpecWidth := logSdSpecWidth]

  rm(dt.sdSpecWidth)
  gc()


  # Tailing factor comparisons ----------------------------

  # create table for max Tfact
  dt.maxTfact <- dt.PkTfact.long[,max(pkTfact, na.rm = TRUE),
                                 keyby = .(featID)]
  colnames(dt.maxTfact)[2] <- "maxTfact"

  # get log maxTfact
  dt.maxTfact[, logMaxTfact := log(maxTfact)]

  # create table for median Tfact
  dt.medTfact <- dt.PkTfact.long[,median(pkTfact, na.rm = TRUE),
                                 keyby = .(featID)]
  colnames(dt.medTfact)[2] <- "medTfact"

  # get log medTfact
  dt.medTfact[, logMedTfact := log(medTfact)]

  # create table for minimum Tfact
  dt.minTfact <- dt.PkTfact.long[,min(log(pkTfact), na.rm = TRUE),
                                 keyby = .(featID)]
  colnames(dt.minTfact)[2] <- "logMinTfact"

  # create table for mean Tfact
  dt.avgTfact <- dt.PkTfact.long[,mean(pkTfact, na.rm = TRUE),
                                 keyby = .(featID)]
  colnames(dt.avgTfact)[2] <- "avgTfact"

  # create table for standard deviation from Tfact
  dt.sdTfact <- dt.PkTfact.long[,sd(pkTfact, na.rm = TRUE),
                                keyby = .(featID)]
  colnames(dt.sdTfact)[2] <- "sdTfact"

  # create log sd Tfact column
  dt.sdTfact[, logSdTfact := log(sdTfact)]

  # add mean to standard deviation table
  dt.sdTfact[dt.avgTfact, avgTfact := avgTfact]

  # create Coefficient of variation column
  dt.sdTfact[, CvTfact := sdTfact/avgTfact]

  # create log cv column
  dt.sdTfact[, logCvTfact := log(sdTfact/avgTfact)]

  # add maxTfact to dt.sdTfact
  dt.sdTfact[dt.maxTfact, maxTfact := maxTfact]

  # add column for sdTfact / maxTfact
  dt.sdTfact[, custVarMeasTfact := sdTfact / maxTfact]

  # add to dt.featSummary
  setkey(dt.featSummary, featID)
  dt.featSummary[dt.maxTfact, logMaxTfact := logMaxTfact]
  dt.featSummary[dt.medTfact, logMedTfact := logMedTfact]
  dt.featSummary[dt.minTfact, logMinTfact := logMinTfact]
  dt.featSummary[dt.sdTfact, sdTfact := sdTfact]
  dt.featSummary[dt.sdTfact, logSdTfact := logSdTfact]
  dt.featSummary[dt.sdTfact, CvTfact := CvTfact]
  dt.featSummary[dt.sdTfact, logCvTfact := logCvTfact]
  dt.featSummary[dt.sdTfact, custVarMeasTfact := custVarMeasTfact]

  rm(dt.maxTfact)
  rm(dt.minTfact)
  rm(dt.sdTfact)
  gc()

  # Asymmetry factor comparisons ----------------------------

  # create table for max Afact
  dt.maxAfact <- dt.PkAfact.long[,max(pkAfact, na.rm = TRUE),
                                 keyby = .(featID)]
  colnames(dt.maxAfact)[2] <- "maxAfact"

  # get log maxAfact
  dt.maxAfact[, logMaxAfact := log(maxAfact)]

  # create table for median Afact
  dt.medAfact <- dt.PkAfact.long[,median(pkAfact, na.rm = TRUE),
                                 keyby = .(featID)]
  colnames(dt.medAfact)[2] <- "medAfact"

  # get log medAfact
  dt.medAfact[, logMedAfact := log(medAfact)]

  # create column for factor by which longer tail is longer than
  # shorter tail in dt.pkAfact.long
  dt.PkAfact.long[, AsymLngFact := ifelse(pkAfact < 1, -1/pkAfact , pkAfact)]

  # create table for max absolute AsymLngFact
  dt.maxAbsAsymLngFact <- dt.PkAfact.long[,max(abs(AsymLngFact), na.rm = TRUE),
                                          keyby = .(featID)]
  colnames(dt.maxAbsAsymLngFact)[2] <- "maxAbsAsymLngFact"

  # add column for logMaxAbsAsymLngFact
  dt.maxAbsAsymLngFact[, logMaxAbsAsymLngFact := log(maxAbsAsymLngFact)]

  # create table for median absolute AsymLngFact
  dt.medLogAbsAsymLngFact <- dt.PkAfact.long[,median(log(abs(AsymLngFact)), na.rm = TRUE),
                                             keyby = .(featID)]
  colnames(dt.medLogAbsAsymLngFact)[2] <- "medLogAbsAsymLngFact"

  # create table for minimum absolute AsymLngFact
  dt.minAbsAsymLngFact <- dt.PkAfact.long[,min(abs(AsymLngFact), na.rm = TRUE),
                                          keyby = .(featID)]
  colnames(dt.minAbsAsymLngFact)[2] <- "minAbsAsymLngFact"

  # create column for log minAbsAsymLngFact
  dt.minAbsAsymLngFact[,logMinAbsAsymLngFact := log(minAbsAsymLngFact)]

  # create table for mean AsymLngFact
  dt.avgAbsAsymLngFact <- dt.PkAfact.long[,mean(abs(AsymLngFact), na.rm = TRUE),
                                          keyby = .(featID)]
  colnames(dt.avgAbsAsymLngFact)[2] <- "avgAbsAsymLngFact"

  # create table for standard deviation of AsymLngFact
  dt.sdAsymLngFact <- dt.PkAfact.long[,sd(AsymLngFact, na.rm = TRUE),
                                      keyby = .(featID)]
  colnames(dt.sdAsymLngFact)[2] <- "sdAsymLngFact"

  # add log standard deviation
  dt.sdAsymLngFact[, logSdAsymLngFact := log(sdAsymLngFact)]

  # add mean to dt.sdAsymLngFact
  dt.sdAsymLngFact[dt.avgAbsAsymLngFact, avgAbsAsymLngFact := avgAbsAsymLngFact]

  # create coefficient of variation column
  dt.sdAsymLngFact[, CvAsymLngFact := sdAsymLngFact/avgAbsAsymLngFact]

  # create log cv column
  dt.sdAsymLngFact[, logCvAsymLngFact := log(CvAsymLngFact)]

  # add maxAbsAsymLngFact to dt.sdAsymLngFact
  dt.sdAsymLngFact[dt.maxAbsAsymLngFact,
                   maxAbsAsymLngFact := maxAbsAsymLngFact]

  # create a column for sdAsymLngFact / maxAbsAsymLngFact
  dt.sdAsymLngFact[, custVarMeasAsymLngFact :=
                     sdAsymLngFact / maxAbsAsymLngFact]

  # add to dt.featSummary
  setkey(dt.featSummary, featID)
  dt.featSummary[dt.maxAfact, logMaxAfact := logMaxAfact]
  dt.featSummary[dt.medAfact, logMedAfact := logMedAfact]
  dt.featSummary[dt.maxAbsAsymLngFact, logMaxAbsAsymLngFact := logMaxAbsAsymLngFact]
  dt.featSummary[dt.medLogAbsAsymLngFact, medLogAbsAsymLngFact := medLogAbsAsymLngFact]
  dt.featSummary[dt.minAbsAsymLngFact, logMinAbsAsymLngFact := logMinAbsAsymLngFact]
  dt.featSummary[dt.sdAsymLngFact, sdAsymLngFact := sdAsymLngFact]
  dt.featSummary[dt.sdAsymLngFact, logSdAsymLngFact := logSdAsymLngFact]
  dt.featSummary[dt.sdAsymLngFact, CvAsymLngFact := CvAsymLngFact]
  dt.featSummary[dt.sdAsymLngFact, logCvAsymLngFact := logCvAsymLngFact]
  dt.featSummary[dt.sdAsymLngFact, custVarMeasAsymLngFact := custVarMeasAsymLngFact]

  rm(dt.maxAfact)
  rm(dt.maxAbsAsymLngFact)
  rm(dt.medLogAbsAsymLngFact)
  rm(dt.minAbsAsymLngFact)
  rm(dt.avgAbsAsymLngFact)
  rm(dt.sdAsymLngFact)
  gc()

  # Tfact / Afact -----------------------------------------------------------

  # add Afact to dt.PkTfact.long
  setkeyv(dt.PkTfact.long, c("featID", "injection"))
  setkeyv(dt.PkAfact.long, c("featID", "injection"))
  dt.PkTfact.long[dt.PkAfact.long, pkAfact := pkAfact]

  # create a column for Tfact / Afact
  dt.PkTfact.long[, pkT_Afact := pkAfact / pkTfact]

  # create table for max T_Afact
  dt.maxT_Afact <- dt.PkTfact.long[,max(pkT_Afact, na.rm = TRUE),
                                   keyby = .(featID)]
  colnames(dt.maxT_Afact)[2] <- "maxT_Afact"

  # get log maxT_Afact
  dt.maxT_Afact[, logMaxT_Afact := log(maxT_Afact)]

  # create table for median T_Afact
  dt.medT_Afact <- dt.PkTfact.long[,median(pkT_Afact, na.rm = TRUE),
                                   keyby = .(featID)]
  colnames(dt.medT_Afact)[2] <- "medT_Afact"

  # get log medT_Afact
  dt.medT_Afact[, logMedT_Afact := log(medT_Afact)]

  # create table for minimum T_Afact
  dt.minT_Afact <- dt.PkTfact.long[,min(log(pkT_Afact), na.rm = TRUE),
                                   keyby = .(featID)]
  colnames(dt.minT_Afact)[2] <- "logMinT_Afact"

  # create table for mean T_Afact
  dt.avgT_Afact <- dt.PkTfact.long[,mean(pkT_Afact, na.rm = TRUE),
                                   keyby = .(featID)]
  colnames(dt.avgT_Afact)[2] <- "avgT_Afact"

  # create table for standard deviation from T_Afact
  dt.sdT_Afact <- dt.PkTfact.long[,sd(pkT_Afact, na.rm = TRUE),
                                  keyby = .(featID)]
  colnames(dt.sdT_Afact)[2] <- "sdT_Afact"

  # create log sd T_Afact column
  dt.sdT_Afact[, logSdT_Afact := log(sdT_Afact)]

  # add mean to standard deviation table
  dt.sdT_Afact[dt.avgT_Afact, avgT_Afact := avgT_Afact]

  # create Coefficient of variation column
  dt.sdT_Afact[, CvT_Afact := sdT_Afact/avgT_Afact]

  # create log cv column
  dt.sdT_Afact[, logCvT_Afact := log(sdT_Afact/avgT_Afact)]

  # add maxT_Afact to dt.sdT_Afact
  dt.sdT_Afact[dt.maxT_Afact, maxT_Afact := maxT_Afact]

  # add column for sdT_Afact / maxT_Afact
  dt.sdT_Afact[, custVarMeasT_Afact := sdT_Afact / maxT_Afact]

  # add to dt.featSummary
  setkey(dt.featSummary, featID)
  dt.featSummary[dt.maxT_Afact, logMaxT_Afact := logMaxT_Afact]
  dt.featSummary[dt.medT_Afact, logMedT_Afact := logMedT_Afact]
  dt.featSummary[dt.minT_Afact, logMinT_Afact := logMinT_Afact]
  dt.featSummary[dt.sdT_Afact, sdT_Afact := sdT_Afact]
  dt.featSummary[dt.sdT_Afact, logSdT_Afact := logSdT_Afact]
  dt.featSummary[dt.sdT_Afact, CvT_Afact := CvT_Afact]
  dt.featSummary[dt.sdT_Afact, logCvT_Afact := logCvT_Afact]
  dt.featSummary[dt.sdT_Afact, custVarMeasT_Afact := custVarMeasT_Afact]


  # Other combined factors ---------------------------------------------------

  # get minFWHM * medTFact * MedAFact
  # based on identified trend in low probability peaks that looked good

  # add medTFact and MedAFact to dt.minRelFWHM
  dt.minRelFWHM[dt.medTfact, medTfact := medTfact]
  dt.minRelFWHM[dt.medAfact, medAfact := medAfact]

  # create column for products
  dt.minRelFWHM[, combFact1 := minRelFWHM * medTfact * medAfact] # more correlated to medAfact than next
  # dt.minRelFWHM[, combFact2 := minRelFWHM * medAfact] # highly correlated with combFact1 lower effect size
  # dt.minRelFWHM[, combFact3 := minRelFWHM * medTfact] # highly correlated with combFact2 lower effect size than 1
  # dt.minRelFWHM[, combFact4 := medTfact * medAfact] # very highly correlated with both medTFact and medAFact

  # get log values
  dt.minRelFWHM[, logCombFact1 := log(combFact1)]
  # dt.minRelFWHM[, logCombFact2 := log(combFact2)]
  # dt.minRelFWHM[, logCombFact3 := log(combFact3)]
  # dt.minRelFWHM[, logCombFact4 := log(combFact4)]

  # variable: custVarMeasRelArea / custVarMeasDur

  # add custVarMeasDur to dt.sdRelArea
  dt.sdRelArea[dt.sdDur, custVarMeasDur := custVarMeasDur]
  dt.sdRelArea[, combVarMeas1 := custVarMeasRelArea / custVarMeasDur]
  dt.sdRelArea[, combVarMeas2 := log(custVarMeasRelArea / custVarMeasDur)]

  # add to dt.featSummary
  setkey(dt.featSummary, featID)
  dt.featSummary[dt.minRelFWHM, logCombFact1 := logCombFact1]
  # dt.featSummary[dt.minRelFWHM, logCombFact2 := logCombFact2]
  # dt.featSummary[dt.minRelFWHM, logCombFact3 := logCombFact3]
  # dt.featSummary[dt.minRelFWHM, logCombFact4 := logCombFact4]
  dt.featSummary[dt.sdRelArea, combVarMeas1 := combVarMeas1]
  dt.featSummary[dt.sdRelArea, combVarMeas2 := combVarMeas2]


  rm(dt.minRelFWHM)
  rm(dt.medTfact)
  rm(dt.medAfact)
  gc()

  print("model variables generated")
  print(runTime(startTime))

  return(dt.featSummary)

  # remove long form tables?

} # end function
