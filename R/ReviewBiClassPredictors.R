#' Generate Tables and Plots Examining Individual Predictors
#'
#' !! This function is in development
#' Creates a new folder in the output directory and generates useful tables
#' and plots for reviewing the predictive value of individual variables
#'
#' @param dt.featSummary.train
#'
#' @return !!!not sure yet
#' @export
#'
#' @examples TODO

ReviewBiClassPredictors <- function(dt.featSummary.train, vect.predCols,
                                    classCol, outputDirPath) {

  # create output directory --------------------------------------

  OutputFolderName <- "ReviewIndividualPredictors"

  outputFolderPath <- CreateOutputDir(outputDirPath, OutputFolderName)

  # examine input table ----------------------------------------------

  # number of observations
  numObs <- nrow(dt.featSummary.train)

  # number of predictors
  numPotPred <- length(vect.predCols)


  # variable correlations (moved for data combination) -------------------------

  # create variable correlation matrix
  dt.corMat <- cor(dt.featSummary.train[,vect.predCols, with = FALSE])
  dt.corMat <- data.table(dt.corMat, keep.rownames = TRUE)

  # convert to long form
  dt.cor.long <- melt(dt.corMat, id.vars = "rn",
                      measure.vars = c(2:ncol(dt.corMat)),
                      variable.name = "corVar",
                      value.name = "corr",
                      na.rm = TRUE)

  # remove correlation of variables with themselves
  dt.cor.long <- subset(dt.cor.long,
                        rn != corVar)

  # add column for absolute value of correlation
  dt.cor.long[, absCorr := abs(corr)]

  # export using defined function
  WrtTable(dt.cor.long, outputFolderPath, "VariableCorrelations")

  print("export variable correlations")

  # review individual variables ---------------------------------

  # create tracking vectors
  vect.varPVal <- numeric(length = numPotPred)
  vect.relAbsEff <- numeric(length = numPotPred)
  vect.mostCorVar <- character(length = numPotPred)
  vect.mostCor <- numeric(length = numPotPred)
  vect.varSd <- numeric(length = numPotPred)
  vect.VarAbsEffSz <- numeric(length = numPotPred)
  vect.avgAbsEr <- numeric(length = numPotPred)
  vect.missIDrate <- numeric(length = numPotPred)

  # create list of empty lists for each page (1 page per predictor)
  plotPage <- rep(list(list()), numPotPred)

  # get name of class column
  nameClassCol <- colnames(dt.featSummary.train)[classCol]

  # formula start text
  lgRegFmla <- paste0(nameClassCol, "  ~ ")

  # store current warning setting
  setWarn <- options()[["warn"]]

  # set warn option to 0 setting
  options(warn = 0)

  # for each variable column, get a linear model p-value
  # and generate a violin plot
  for (varNum in 1:numPotPred) {

    # variable column
    lpvar.colNum <- vect.predCols[varNum]

    # variable column name
    lpvar.varName <- colnames(dt.featSummary.train)[lpvar.colNum]

    # extract variable column
    lpvar.varCol <- subset(dt.featSummary.train,
                           select = c(lpvar.colNum))

    # variable min
    lpvar.varMin <- min(lpvar.varCol)

    # variable max
    lpvar.varMax <- max(lpvar.varCol)

    # create table from variable and class columns
    lptbl.varCol <- subset(dt.featSummary.train,
                           select = c(lpvar.colNum,
                                      (ncol(dt.featSummary.train) - 1),
                                      ncol(dt.featSummary.train)))

    # rename variable column for easy use
    colnames(lptbl.varCol)[1] <- "var"

    # variable standard deviation
    lpvar.sd <- sd(lptbl.varCol$var, na.rm = TRUE)

    # table of means by group
    lptbl.grpMeans <- lptbl.varCol[, mean(var, na.rm = TRUE),
                                   keyby = .(class)]
    colnames(lptbl.grpMeans)[2] <- "avgVar"

    # absolute effect size
    lpvar.absEffsz <- abs(lptbl.grpMeans$avgVar[1] - lptbl.grpMeans$avgVar[2])

    # relative absolute effect size
    lpvar.relAbsEffsz <- lpvar.absEffsz / lpvar.sd

    # subset correlation table
    lptbl.cor <- subset(dt.cor.long,
                        rn == lpvar.varName)

    # subset to max absolute correlation
    lptbl.cor <- subset(lptbl.cor,
                        absCorr == max(absCorr))

    # get most correlated variable
    lpvar.mostCorVar <- as.character(lptbl.cor$corVar)

    # get correlation for most correlated
    lpvar.mostCor <- round(lptbl.cor$corr,4)
    lpvar.mostCorPrnt <- round(lpvar.mostCor,4)

    # text for formula
    lpvar.formText <- paste0(lgRegFmla,lpvar.varName)

    # create formula object
    lpvar.fmla <- as.formula(lpvar.formText)

    # create logistic regression model
    lpvar.lm <- glm(lpvar.fmla,
                    family = binomial,
                    data = dt.featSummary.train)

    # use predict to get a vector of the predicted probablilities of good
    lpvect.predCheck <- predict(lpvar.lm,
                                newdata = dt.featSummary.train,
                                type = "response")

    # add predicted probability
    lptbl.varCol[, predProb := lpvect.predCheck]

    # create column for error
    lptbl.varCol[, error := numClass - predProb]

    # create column for absolute error
    lptbl.varCol[, absError := abs(error)]

    # create column for misidentification (absolute error > 0.5)
    lptbl.varCol[, missID := ifelse(absError > 0.5, 1, 0)]

    # find average absolute error
    lpvar.avgAbsEr <- mean(lptbl.varCol$absError)

    # calculate missidentification rate
    lpvar.missIDrate <- mean(lptbl.varCol$missID)

    # get p-value
    lpvar.pval <- coef(summary(lpvar.lm))[2,4]
    # format for printing
    lpvar.pvalprnt <- sprintf("%1.2E",lpvar.pval)

    # add to tracking vectors
    vect.mostCorVar[varNum] <- lpvar.mostCorVar
    vect.mostCor[varNum] <- lpvar.mostCor
    vect.varSd[varNum] <- lpvar.sd
    vect.VarAbsEffSz[varNum] <- lpvar.absEffsz
    vect.relAbsEff[varNum] <- lpvar.relAbsEffsz
    vect.varPVal[varNum] <- lpvar.pval
    vect.avgAbsEr[varNum] <- lpvar.avgAbsEr
    vect.missIDrate[varNum] <- lpvar.missIDrate

    # create plot title
    lpvar.infoTitle <- paste0(" P-value = ", lpvar.pvalprnt,
                              "\n Rel Abs Effect = ", lpvar.relAbsEffsz,
                              "\n Most correlated: ", lpvar.mostCorVar,
                              " ", lpvar.mostCorPrnt,
                              "\n AvgAbsError: ", lpvar.avgAbsEr)

    # violin plots
    plotPage[[varNum]][[1]] <-
      ggplot() +
      geom_violin(data = dt.featSummary.train,
                  aes_string(x = "class", y = lpvar.varName),
                  fill = NA, size = .75, draw_quantiles = c(0.5)) +
      # geom_point(data = dt.featSummary.train, aes_string(x = "class", y = lpvar.varName), shape = ".") +
      ggtitle(lpvar.varName) +
      ylab(lpvar.varName) +
      xlab('Classification') +
      ggThemeLcmsMetab()

    # Logistic Regression Plots

    # prediction values for creating curve
    lptbl.plotPred <- data.table(Xpred = seq(from = lpvar.varMin,
                                             to = lpvar.varMax,
                                             length.out = 200))
    colnames(lptbl.plotPred)[1] <- lpvar.varName

    # prediction values vectors
    lpvect.plotYpred <- predict(lpvar.lm,
                                newdata = lptbl.plotPred,
                                type = "response")

    # add to plot data table
    lptbl.plotPred[, probPred := lpvect.plotYpred]

    # logistic regression plots
    plotPage[[varNum]][[2]] <-
      ggplot() +
      # geom_point(data = dt.featSummary.train,
      # aes_string(x = lpvar.varName, y = "numClass"), shape = ".") +
      geom_line(data = lptbl.plotPred,
                aes_string(x = lpvar.varName, y = "probPred")) +
      geom_hline(yintercept = 1) +
      geom_hline(yintercept = 0) +
      ggtitle(lpvar.infoTitle) +
      ylab("predicted probability") +
      xlab(lpvar.varName) +
      ggThemeLcmsMetab()

  } # end variable examination for loop with plots

  # return warning option
  options(warn = setWarn)

  # create data table for variable summary
  dt.varSummary <- data.table(var = colnames(dt.featSummary.train)[4:(3 + numPotPred)],
                              mostCorVar = vect.mostCorVar,
                              mostCor = vect.mostCor,
                              absMostCor = abs(vect.mostCor),
                              StdDev = vect.varSd,
                              absEffSz =  vect.VarAbsEffSz,
                              relAbsEff = vect.relAbsEff,
                              pVal_1Regr = vect.varPVal,
                              avgAbsEr = vect.avgAbsEr,
                              missIDrate = vect.missIDrate)

  logText <-
    UpdateLogText(logText,"individual variable review completed")
  UpdateLogFile(logFilePath, logText)

  # PDF individual variable plots --------------------------------------------------

  PdfFilePath <- paste0(outputFolderPath,"/",
                        "IndividualVariablePlots",
                        startTimeStamp,".pdf")

  pdf(PdfFilePath, width = 10, height = 7.5, paper = "USr", onefile = TRUE)
  for (i in seq(length(plotPage))) {
    grid.arrange(grobs = plotPage[[i]], ncol = 2)

    percComplete <- i/length(plotPage) * 100

    print(paste0("individual variable review pdf completion ",
                 format(percComplete, nsmall = 1),"%"))
  }
  dev.off()

  logText <-
    UpdateLogText(logText,"individual variable review PDF file generated")
  UpdateLogFile(logFilePath, logText)

  # plot avgAbsError & missIDrate Individual variables ---------------------------------------------------------

  # add mean avgAbsErRank  as additional column
  dt.varSummary[, avgAbsErRank := frank(-avgAbsEr)]

  # pass the ranking on to the variable names
  # dt.varSummary[, var := reorder(var,avgAbsErRank)]

  # sort by rank for plot labelling
  setkey(dt.varSummary, avgAbsErRank)

  # convert avgAbsErRank to factor
  dt.varSummary$avgAbsErRank <- as.factor(dt.varSummary$avgAbsErRank)

  # save filenames in order as vector
  vect.varNames <- dt.varSummary$var

  # create empty plotlist
  plotVarModEr <- list()
  plotVarModEr[[1]] <- list()



  # create plot object average absolute error
  plotVarModEr[[1]][[1]] <-
    ggplot() +
    geom_point(data = dt.varSummary,
               aes(x = avgAbsErRank, y = round(avgAbsEr,3))) +
    ggtitle("Single Variable Model Error") +
    ylab("average absolute error") +
    xlab("prediction variable") +
    scale_x_discrete(labels = vect.varNames) +
    ggThemeLcmsMetab()

  plotVarModEr[[2]] <- list()

  # create plot object average absolute error
  plotVarModEr[[2]][[1]] <-
    ggplot() +
    geom_point(data = dt.varSummary,
               aes(x = avgAbsErRank, y = round(missIDrate,3))) +
    ggtitle("Single Variable Model Misidentification Rate") +
    ylab("misidentification rate") +
    xlab("prediction variable") +
    scale_x_discrete(labels = vect.varNames) +
    ggThemeLcmsMetab()

  # create pdf
  PdfPlotsSingleFile(plotVarModEr, "landscape", outputFolderPath, "PlotIndividualVarAvgAbsErMissID")

  logText <-
    UpdateLogText(logText,"plot complete: individual avg abs error and missidentification rate")
  UpdateLogFile(logFilePath, logText)


} # end function
