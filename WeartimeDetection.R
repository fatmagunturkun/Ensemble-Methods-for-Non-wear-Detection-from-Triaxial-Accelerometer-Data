# This function runs weartime detection using the Choi algorithm
# on the specified accelerometer data
weartimeChoi <- function(acc, usefixed,frame = 90, streamFrame = 30)
{
    # if it is required to use fixedChoi, then load source, otherwise load package
    # and load code with the fix
    if (usefixed) {
        source('choiFixed.R')
    } else {
        library("PhysicalActivity")
    }
    
    # check appropriate fields
    checkActigraphField(acc, 'Date')
    checkActigraphField(acc, 'Time')
    checkActigraphField(acc, 'Axis1')
    
    # generate time stamp
    acc$date_form <- mdy(acc$Date)
    acc$time_form <- hms(acc$Time)

    acc$TimeStamp <- paste0(year(acc$date_form), "-",
    str_pad(month(acc$date_form), width = 2, side = "left", pad = "0"), "-",
    str_pad(day(acc$date_form), width = 2, side = "left", pad = "0"), " ",
    str_pad(hour(acc$time_form), width = 2, side = "left", pad = "0"), ":",
    str_pad(minute(acc$time_form), width = 2, side = "left", pad = "0"), ":",
    str_pad(second(acc$time_form), width = 2, side = "left", pad = "0"))

    # if fixed Choi was used, functions, otherwise, unlink the package
    if (usefixed) {
        
        # extract table to use for weartime detection
        pr_acc <- acc[ ,c("TimeStamp", "Axis1")]
        colnames(pr_acc)[colnames(pr_acc) == "Axis1"] <- "counts"
        # apply Choi
        wearmark2X <- wearingMarking(pr_acc, frame = frame, perMinuteCts = 60, streamFrame = streamFrame)
        
        pr_acc <- acc[ ,c("TimeStamp", "Axis2")]
        colnames(pr_acc)[colnames(pr_acc) == "Axis2"] <- "counts"
        wearmark2Y <- wearingMarking(pr_acc, frame = frame, perMinuteCts = 60, streamFrame = streamFrame)
        
        pr_acc <- acc[ ,c("TimeStamp", "Axis3")]
        colnames(pr_acc)[colnames(pr_acc) == "Axis3"] <- "counts"
        wearmark2Z <- wearingMarking(pr_acc, frame = frame, perMinuteCts = 60, streamFrame = streamFrame)
        
        pr_acc <- acc[ ,c("TimeStamp", "Vector.Magnitude")]
        colnames(pr_acc)[colnames(pr_acc) == "Vector.Magnitude"] <- "counts"
        wearmark2VM <- wearingMarking(pr_acc, frame = frame, perMinuteCts = 60, streamFrame = streamFrame)
        rm(list='wearingMarking', envir=sys.frame(-2))
        rm(list='dataCollapser', envir=sys.frame(-2))
        rm(list='markingTime', envir=sys.frame(-2))
    } else {
      # apply Choi
        wearmark2X <- wearingMarking(dataset = acc, TS = "TimeStamp", cts = "Axis1", frame = frame, perMinuteCts = 60, streamFrame = streamFrame)
        wearmark2Y <- wearingMarking(dataset = acc, TS = "TimeStamp", cts = "Axis2", frame = frame, perMinuteCts = 60, streamFrame = streamFrame)
        wearmark2Z <- wearingMarking(dataset = acc, TS = "TimeStamp", cts = "Axis3", frame = frame, perMinuteCts = 60, streamFrame = streamFrame)
        wearmark2VM <- wearingMarking(dataset = acc, TS = "TimeStamp", cts = "Vector.Magnitude", frame = frame, perMinuteCts = 60, streamFrame = streamFrame)
        detach("package:PhysicalActivity", unload=TRUE)
    }
    # check that the results are aligned with the inputs in terms of timestamps
    if(sum(wearmark2X$TimeStamp != acc$TimeStamp)>0) {
      acc$TimeStampW <- wearmark2X$TimeStamp
      warning('Algorithm produced mismatched timestamps X. Please check TimeStampW field.')
    }
    # check that the results are aligned with the inputs in terms of timestamps
    if(sum(wearmark2Y$TimeStamp != acc$TimeStamp)>0) {
      warning('Algorithm produced mismatched timestamps Y. Please check TimeStampW field.')
    }
    # check that the results are aligned with the inputs in terms of timestamps
    if(sum(wearmark2Z$TimeStamp != acc$TimeStamp)>0) {
      warning('Algorithm produced mismatched timestamps Z. Please check TimeStampW field.')
    }
    # check that the results are aligned with the inputs in terms of timestamps
    if(sum(wearmark2VM$TimeStamp != acc$TimeStamp)>0) {
      warning('Algorithm produced mismatched timestamps VM. Please check TimeStampW field.')
    }
    # record the wear marking output
    acc$weartimeChoiX <- as.numeric(wearmark2X$wearing=='w')
    acc$weartimeChoiY <- as.numeric(wearmark2Y$wearing=='w')
    acc$weartimeChoiZ <- as.numeric(wearmark2Z$wearing=='w')
    acc$weartimeChoiVM <- as.numeric(wearmark2VM$wearing=='w')
    
    acc$weekday <- wearmark2X$weekday
    acc$days <- wearmark2X$days
    return (acc)
}

# The function checks if the specified field is present in the acc data frame
checkActigraphField <- function(acc, fieldname)
{
    if (!with(acc, exists(fieldname))) {
        stop(paste('Actigraph data is missing the', fieldname, 'field. Please check that the file is read correctly.', sep = ' '))
    }
}

# This function runs weartime detection using the NHANES algorithm
# We apply it to the 3 axes and the step counts
weartimeNHANES <- function(acc,window, tol, tol_upper, nci, days_distinct)
{
  # set up workspace
  library("accelerometry")
  
  # set number of seconds to be aggregated for nhanes analysis (60)
  nsec = 60
  
  # convert the accelerometry data to 1 minute intervals
  nhanes = data.frame(
    Y = convertNHANES(acc, 'Axis1', nsec, 1),
    X = convertNHANES(acc, 'Axis2', nsec, 1),
    Z = convertNHANES(acc, 'Axis3', nsec, 1),
    VM = convertNHANES(acc, 'Vector.Magnitude', nsec, 1))
  
  # detect the weartime using the nhanes algorithm
  # alternative: tol = 0, tol.upper = 99, nci = FALSE
  # alternative: tol = 1, tol.upper=99, nci = TRUE
  
  nhanes$weartimeX <- weartime(nhanes$X, window = window, tol = tol, tol_upper=tol_upper, nci = nci, days_distinct = days_distinct)
  nhanes$weartimeY <- weartime(nhanes$Y, window = window, tol = tol, tol_upper=tol_upper, nci = nci, days_distinct = days_distinct)
  nhanes$weartimeZ <- weartime(nhanes$Z, window = window, tol = tol, tol_upper=tol_upper, nci = nci, days_distinct = days_distinct)
  nhanes$weartimeVM <- weartime(nhanes$VM, window = window, tol = tol, tol_upper=tol_upper, nci = nci, days_distinct = days_distinct)
  
  # revert to 1 second intervals by duplicating weartime
  nhanes <- nhanes[rep(seq_len(nrow(nhanes)), each=nsec),]
  nhanes <- nhanes[1:nrow(acc),]
  
  # copy the weartime to the acc data frame
  acc$weartimeNhanesX <- nhanes$weartimeX
  acc$weartimeNhanesY <- nhanes$weartimeY
  acc$weartimeNhanesZ <- nhanes$weartimeZ
  acc$weartimeNhanesVM <- nhanes$weartimeVM
  
  # clean up workspace
  detach("package:accelerometry", unload=TRUE)
  
  return (acc)
}

# Convert 1 second intervals to larger intervals, of size specified by numsec,
# and makes a matrix with numcol elements per row
convertNHANES <- function(acc, fieldname, numsec, numcol)
{
  counts1sec <- acc[,fieldname]
  d = length(counts1sec)
  d1 = numsec # should be 60 for aggregation to 1 minute
  d2 = numcol # should be 60 if d1 is 60, to have 1 hour represented per row, where each col represents 1 minute aggregates
  d3 = ceiling(length(counts1sec)/(d1*d2))
  if (length(counts1sec) > d1*d2*d3) {
    counts1sec <- counts1sec[1:(d1*d2*d3),]
  }
  if (length(counts1sec) < d1*d2*d3) {
    counts1sec[(d+1):(d1*d2*d3)] = 0
  }
  dim(counts1sec) <- c(d1,d2,d3)
  # If M is a 3D matrix and d1=dim(M)[1], d2=dim(M)[2] then
  #    M[i,j,k]==M[ i+ (j-1)*d1 + (k-1)*d1*d2 ]
  nhanes <- colSums(counts1sec, dims = 1)
  nhanes <- t(nhanes)
  return (nhanes)
}

# make vector v of the specified length n by truncating or padding with 1s.
shapevector <- function(v, n) {
    if (length(v)<n) {
        p = rep.int(c(1), n-length(v))
        v <- append(v,p)
    }
    if (length(v)>n) {
        v <- v[1:n]
    }
    return (v)
}

# Compute weartime through HMM 
weartimeHMM <- function(acc) {
  g <- as.integer(3)  # number of observations
  w <- as.integer(1200)  # size of window for average computation
  
  # Load reticulate library in R
  library(reticulate)
  
  # Use Python within R
  numpy <- import("numpy")
  pandas <- import("pandas")
  hmmlearn <- import("hmmlearn.hmm")
  
  # compute the mean absolute value of the acc data over each time window
  macc <- data.frame(matrix(ncol = g, nrow = floor(nrow(acc)/w)))
  for (i in 0:floor(nrow(acc)/w)) {
    range <- ((w*i)+1):min(w*(i+1), nrow(acc))
    X <- acc[range, c("Axis1", "Axis2", "Axis3")]
    X[is.na(X)] <- 0
    macc[i+1,] <- colMeans(abs(X))
  }
  
  # Number of states
  n_states_py <- as.integer(4)
  
  # Function to fit model and return weartime
  fit_model <- function(model_type, scale_flag) {
    model <- hmmlearn$GaussianHMM(n_components = n_states_py, covariance_type = model_type, n_iter = as.integer(1000))#init_params = 'random'
    
    if (scale_flag) {
      # Create pandas DataFrame
      macc_py <- r_to_py(scale(macc), convert = TRUE)
      model$fit(macc_py)
      hidden_states <- model$predict(macc_py)
    } else {
      macc_py <- r_to_py(macc, convert = TRUE)
      model$fit(macc_py)
      hidden_states <- model$predict(macc_py)
    }
    
    # obtain the average magnitude over each cluster
    magnitude <- rep(0.0, n_states_py)
    for (k in 0:(n_states_py-1)) {
      Xk <- macc[hidden_states == k,]
      magnitude[k+1] <- mean(colMeans(abs(Xk)))
    }
    knw <- which.min(magnitude)
    weartime <- rep(1, nrow(macc))
    weartime[hidden_states == (knw-1)] <- 0
    return(weartime)
  }
  # Fit model
  weartime <- fit_model('spherical',FALSE)
  
  acc$weartimeHMM <- 0
  # for each window, assign all the points in the window the wear classification corresponding to the window
  for (i in 0:floor(nrow(acc)/w)) {
    range <- ((w*i)+1):min(w*(i+1), nrow(acc))
    acc$weartimeHMM[range] <- weartime[i+1]
  }
  return (acc)
}

# This function runs weartime detection on the specified actigraph file
# The file should contain 1 secound counts

runWeartimeDetection <- function(actigraphfile, directory,subid, usefixedchoi,frame, streamFrame, numskiplines,window, tol, tol_upper, nci, days_distinct)
{
  # read actigraph file, skip everything before the header
  print(paste('Processing file ', actigraphfile))
  acc <- read.table(file = actigraphfile, sep = ",", header = T, stringsAsFactors = F, skip = numskiplines)
  # generate and assign timestamps
  ts <- seq(c(ISOdate(2017,1,1)), by = "sec", length.out = nrow(acc))
  acc$Date <- strftime(ts, format="%m/%d/%Y")
  acc$Time <- strftime(ts, format="%H:%M:%S %p")
  
  # rename columns
  acc <- renamecolumn(acc, 'waistwaist_axis1', 'Axis1')
  acc <- renamecolumn(acc, 'waistwaist_axis2', 'Axis2')
  acc <- renamecolumn(acc, 'waistwaist_axis3', 'Axis3')
  acc <- renamecolumn(acc, 'waistwaist_steps', 'Steps')
  acc <- renamecolumn(acc, 'waistwaist_lux', 'Lux')
  acc <- renamecolumn(acc, 'waistwaist_incl_off', 'Inclinometer.Off')
  acc <- renamecolumn(acc, 'waistwaist_incl_stand', 'Inclinometer.Standing')
  acc <- renamecolumn(acc, 'waistwaist_incl_sit', 'Inclinometer.Sitting')
  acc <- renamecolumn(acc, 'waistwaist_incl_lying', 'Inclinometer.Lying')
  acc <- renamecolumn(acc, 'waistwaist_vm', 'Vector.Magnitude')
  
  # run all available weartime detection algorithms
  acc <- weartimeChoi(acc, usefixedchoi,frame=frame, streamFrame=streamFrame)
  acc <- weartimeNHANES(acc,window=window, tol=tol, tol_upper=tol_upper, nci=nci, days_distinct=days_distinct)
  acc <-weartimeHMM(acc)
 
  # write the weartime detection result
  outfile <- paste(directory, paste(toString(subid), '_weartime_1.csv', sep=''), sep="/")
  write.table(acc, file = outfile, sep = ",", row.names = F, col.name = T)
  return (acc)
}
