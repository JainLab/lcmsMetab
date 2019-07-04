#' Generate Model Building Table By Combining Training Data Sets
#'
#' Training data files are provided in a table of file names including
#' the peak attribute tables in the first column and the corresponding tables
#' accepted peaks' mz and RT values in the second column. This function
#' reads the files and generates the table used for model building.
#'
#' @param parentDir path; path to directory with the peak attribute
#' tables of the training data, and the key table
#' @param buildDataKeyFileName file name of the csv table with the peak
#' attribute file name in column 1 and the mz-rt of retained peaks file in
#' column 2 imported above
#' @param outputFolderPath path to directory where build data key will be
#' output
#'
#' @return returns a data table in the format required to build a prediction
#' model
#' @export
#'
#' @examples TODO
#'
CreateBuildData <- function(parentDir, buildDataKeyFileName, outputFolderPath) {

  # import file name table --------------------------

  buildDataKeyFilePath <- paste0(parentDir, "/", buildDataKeyFileName)

  dt.buildTables <- fread(buildDataKeyFilePath, sep = ",",
                          header = TRUE, select = c(1:2))

  # rename columns
  colnames(dt.buildTables) <- c("peakCharFile", "retainedPeakFile")

  # number of build data sets to be importd
  numBuldDataSets <- nrow(dt.buildTables)

  # create vector of randomly generated letter strings to identify
  # each file.
  vect.randStr <- do.call(paste0,
                          replicate(5, sample(LETTERS, numBuldDataSets, TRUE), FALSE))

  # add as new column
  dt.buildTables[, fileID := vect.randStr]

  # export table with IDs added -------------------------------------

  # use defined function
  WrtTable(dt.buildTables, outputFolderPath, "buildDataIDTable")

  # function for reading build data and generating prediction table
  # to be run on each pair of build files
  readAndPredict <- function(pkCharFileBuild, parentDir, dt.buildTables) {

    # import table
    pkCharFileBuildPath <- paste0(parentDir, "/", pkCharFileBuild)

    dt.AllCharTrain <- fread(pkCharFileBuildPath,
                             sep = ",",
                             header = TRUE)

    print(paste0("imported ", pkCharFileBuild))

    # Apply feature summary function
    dt.featSummary.train <- PeakPredictionTable(dt.AllCharTrain)

    print("feature summary table generated")

    # get file id and add to featID ----------------------------

    dt.buildTables.subset <- subset(dt.buildTables,
                                    peakCharFile == pkCharFileBuild)

    fileID <- dt.buildTables.subset$fileID[1]

    # add file id to featID
    dt.featSummary.train[, featID := paste0(featID,"_",fileID)]

    # import classification data -----------------------------------

    # select retained peaks file name from table
    retainedFileTrain <- dt.buildTables.subset$retainedPeakFile[1]

    retainedFileTrainPath <- paste0(parentDir, "/",  retainedFileTrain)

    # import the retained peaks mz-rt table
    dt.retained <- fread(retainedFileTrainPath, sep = ",",
                         select = c(1:2))

    # remove duplicate rows
    dt.retained <- unique(dt.retained)

    # add classification to dt.featSummary.train -----------------

    # create class column of all 1s
    dt.retained[, class := rep.int(TRUE,nrow(dt.retained))]

    # first create featID column in dt.retained
    dt.retained <- AddMzRtFeatIdCol(dt.retained, mzColNum = 1,
                                    rtColNum = 2, placeAfterColNum = 2,
                                    colName = "featID")

    # add file id to featID
    dt.retained[, featID := paste0(featID,"_",fileID)]

    # add class column to dt.featSummary.train
    setkey(dt.retained, featID)
    setkey(dt.featSummary.train, featID)
    dt.featSummary.train[dt.retained, class := class]

    # assign the unretained features FALSE in class column
    dt.featSummary.train[, class := ifelse(is.na(class), FALSE, class)]

    # create numeric class column
    dt.featSummary.train[, numClass := as.numeric(class)]

    return(dt.featSummary.train)

  }

  # lapply function to create feature summary tables for all data
  # as list of tables
  dt.featSummary.build <- lapply(dt.buildTables$peakCharFile,
                                 readAndPredict,
                                 parentDir, dt.buildTables)

  # combine all tables in list
  dt.featSummary.build <- rbindlist(dt.featSummary.build,
                                    use.names = TRUE)

  return(dt.featSummary.build)

}
