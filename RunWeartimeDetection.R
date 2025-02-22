### Main program

rm(list = ls())

library(lubridate)
library(stringr)

source('WeartimeDetection.R')
source('confusion_matrix.R')

renamecolumn <- function(acc, oldname, newname) {
    acc[,newname] <- acc[,oldname]
    acc[,oldname] <- NULL
    return (acc)
}

data_directory <- 'path/to/data'

subjectids = c(list_of_IDs)

directory <-'path/to/outputs'
# Initialize an empty list to store summary data frames

summary_data_list <- list()
start_from_scratch = TRUE

for (subid in subjectids) {
  
  # generate weartime and actigraph file names based on subject id
  actigraphfile <- paste(data_directory,paste('merge1sec_',toString(subid),'_compressed.csv',sep=''), sep="/")
  weartimefile <- paste(directory, paste(toString(subid), '_weartime_1sec.csv', sep=''), sep="/")
  
  # run weartime detection
  if (start_from_scratch) {
    acc <-runWeartimeDetection(actigraphfile, directory,subid, usefixedchoi=TRUE, frame = 90, streamFrame = 30,numskiplines=0, window = 60, tol = 2, tol_upper=100, nci = FALSE, days_distinct = FALSE)
    
    # get ground truth
    acc$weartimeTruth <- 1-as.integer(acc$se_event=='Off')
    
    acc_sum <- acc[c("TimeStamp", "weartimeTruth", "Vector.Magnitude")]
    
    # Calculate time intervals
    acc_sum$interval <- cumsum(c(0, diff(acc_sum$weartimeTruth) != 0))
    
    # Summarize switches and total time
    switches <- table(acc_sum$interval, acc_sum$weartimeTruth)
    total_time <- sum(table(acc_sum$interval))
    
    # Calculate the percentage of 0 values in each interval
    percentage_zeros <- prop.table(table(acc_sum$interval, acc_sum$Vector.Magnitude == 0), margin = 1) * 100
    
    # Combine switches and percentage of 0 values
    summary_data <- data.frame(Switches_0 = switches[,1],Switches_1 = switches[,2], PercentageZeros = percentage_zeros[, 2])
    
    # Extract ID from the file name
    summary_data$id <- subid
    
    # Append the summary_data to the list
    summary_data_list <- c(summary_data_list, list(summary_data))
    
  } else {
    acc <- read.table(file = weartimefile, sep = ",", header = T, stringsAsFactors = F)
    print(paste('Loaded file ', weartimefile))
  }
  
  # write the weartime
  write.table(acc, file = weartimefile, sep = ",", row.names = F, col.name = T)
  
  # retrieve the weartime predictions given by the experts
  wt <- getWeartimePredictions(acc)
  
  # compute the confusion matrix, unrolled to a vector per subject
  confusion <- getConfusionMatrix(wt, wt$weartimeTruth)
  confusion$subid = subid
  
  if (exists("confusion_table")) {
    confusion_table <- rbind(confusion_table, confusion)
    
  } else {
    confusion_table <- confusion 
  }
  
  # plot wear detection graph
  pngfile <- paste(directory, paste(toString(subid), '_wear_detection.png', sep=''), sep="/")
  plotWearGraph(acc, pngfile)
}

# Combine all summary_data frames in the list row-wise
combined_summary_data <- do.call(rbind, summary_data_list)
write.table(combined_summary_data, file = paste(directory,"summary_data.csv", sep="/"), 
            sep = ",", row.names = F, col.name = T)
# write the resulting tables
write.table(confusion_table, file = paste(directory,"confusion_table.csv", sep="/"), 
            sep = ",", row.names = F, col.name = T)

