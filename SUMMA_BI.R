source('WeartimeEnsemble.R')

# Load the installed package
library("summa")
library(tidyr)
library(dplyr)


############################################################
#apply SUMMA and BI 

outputdirectory <- 'path/to/outputs'

# List of models to keep
models_to_keep <- c('weartimeChoiX','weartimeChoiY','weartimeChoiZ','weartimeNhanesX','weartimeNhanesY','weartimeNhanesZ','weartimeHMM')

subjectids = c(list_of_IDs)

for (subid in subjectids) {
  
  
  weartimefile <- paste(outputdirectory, paste(toString(subid), '_weartime_1sec.csv', sep=''), sep="/")
  wt <- read.table(file = weartimefile, sep = ",", header = T)
  
  #apply BI
  input_data <- wt[, models_to_keep]
  
  input_data$ID <- 1:nrow(input_data)
  
  # Convert the data to long format
  input_data <- input_data %>%
    pivot_longer(cols = -ID, names_to = "FUNCTION_ID", values_to = "VALUE") %>%
    mutate(LABEL = "") %>%
    select(ID, LABEL, FUNCTION_ID, VALUE)
  
  # Sort the dataframe by FUNCTION_ID
  input_data <- input_data[order(input_data$FUNCTION_ID), ]
  # Create a mapping of unique function IDs to numeric values
  function_id_mapping <- unique(input_data$FUNCTION_ID)
  names(function_id_mapping) <- function_id_mapping
  
  # Convert FUNCTION_ID values to numeric values
  input_data$FUNCTION_ID <- match(input_data$FUNCTION_ID, names(function_id_mapping))
  
  write.csv(input_data, "input_data.csv", row.names = FALSE, fileEncoding = "UTF-8",quote = FALSE )
  
  system('java -cp "./makina.jar" makina.learn.classification.reflection.Integrator -d input_data.csv -e error_rates.csv -i integrated_data.csv -m BI -o 4000:10:200:-:-:-:-')
  
  # Read the integrated data from the CSV file
  integrated_data <- read.csv("integrated_data.csv")
  
  wt$bi_2 <- integrated_data$VALUE
  wt$bi_2 <- ifelse(wt$bi_2 > 0.5, 1, 0)
  
  # retrieve the weartime predictions given by the experts
  pre <- wt[, models_to_keep]
  pre[pre == 0] <- -1
  
  #apply summa
  summa<-summa(pre,"binary")
  summa@estimated_label[summa@estimated_label == -1] <- 0
  wt$summa_2 <- summa@estimated_label
  
  
  # compute the confusion matrix, unrolled to a vector per subject
  confusion <- getConfusionMatrix(wt[,c("summa_2","bi_2")], wt$weartimeTruth)
  confusion$subid = subid
  write.table(wt, file = paste(outputdirectory, paste(toString(subid), '_weartime_1sec.csv', sep=''), sep="/"), sep = ",", row.names = F, col.name = T)
  
  if (exists("confusion_table")) {
    confusion_table <- rbind(confusion_table, confusion)
    
  } else {
    confusion_table <- confusion
    
  }
  # Remove input data and integrated data files
  file.remove("input_data.csv", "integrated_data.csv","error_rates.csv")
}

# write the resulting tables
write.table(confusion_table, file = paste(outputdirectory,"confusion_table.csv", sep="/"), 
            sep = ",", row.names = F, col.name = T)

##################################################
#plots
library(RColorBrewer)
plotWearGraph <- function(acc, filename)
{
  # Set up the number of plots
  numplots <- 11
  
  # Create custom palettes for each group with 6 colors
  choi_palette <- brewer.pal(6, "BuPu")
  nhanes_palette <- brewer.pal(6, "OrRd")
  zero_counts_palette <- brewer.pal(6, "YlGnBu")
  summa_palette <- c("chartreuse1", "chartreuse2", "chartreuse3", "chartreuse4", "chartreuse3", "chartreuse2")
  bi_palette <- c("cyan1", "cyan2", "cyan3", "cyan4", "cyan3", "cyan2")
  hmm_palette <- brewer.pal(6, "Blues")
  
  # Set up the PNG file
  png(filename, width = 10, height = 12, units="in", res=300)
  par(mfrow=c(numplots,1), mar=c(4, 5, 1, 2))  # Adjust margin
  
  # Plotting Choi data if available
  if (with(acc, exists('weartimeChoiX'))) {
    plot(acc$weartimeChoiX, xlab="", ylab="Choi X", col=choi_palette[2],cex.lab=1.4)
    plot(acc$weartimeChoiY, xlab="", ylab="Choi Y", col=choi_palette[3],cex.lab=1.4)
    plot(acc$weartimeChoiZ, xlab="", ylab="Choi Z", col=choi_palette[4],cex.lab=1.4)
  }
  
  # Plotting NHANES data if available
  if (with(acc, exists('weartimeNhanesX'))) {
    plot(acc$weartimeNhanesX, xlab="", ylab="NHANES X", col=nhanes_palette[2],cex.lab=1.3)
    plot(acc$weartimeNhanesY, xlab="", ylab="NHANES Y", col=nhanes_palette[3],cex.lab=1.3)
    plot(acc$weartimeNhanesZ, xlab="", ylab="NHANES Z", col=nhanes_palette[4],cex.lab=1.3)
  }
  
  # Plotting HMM data if available
  if (with(acc, exists('weartimeHMM_s'))) {
    plot(acc$weartimeHMM_s, xlab="", ylab="HMM", col=hmm_palette,cex.lab=1.4)
  }
  
  # Plotting SUMMA data if available
  if (with(acc, exists('summa'))) {
    plot(acc$summa, xlab="", ylab="SUMMA", col=summa_palette[2],cex.lab=1.4)
    #plot(acc$summa_2, xlab="", ylab="SUMMA2", col=summa_palette[4])
  }
  
  # Plotting BI data if available
  if (with(acc, exists('bi'))) {
    plot(acc$bi, xlab="", ylab="BI", col=bi_palette[2],cex.lab=1.4)
    #plot(acc$bi_2, xlab="", ylab="BI2", col=bi_palette[4])
  }
  
  if (with(acc, exists('weartimeTruth'))) {
    plot(acc$weartimeTruth, xlab="", ylab="Actual", col = "red",cex.lab=1.4)
  }
  
  plot(as.numeric(acc$Axis1), xlab="Seconds", ylab="Counts", col = "black", pch=".",cex.lab=1.4)
  
  title(main="Wear Detection")
  dev.off()
}

for (subid in subjectids) {
  weartimefile <- paste(directory, paste(toString(subid), '_weartime_1sec.csv', sep=''), sep="/")
  acc <- read.table(file = weartimefile, sep = ",", header = T, stringsAsFactors = F, skip = 0)
  #print(acc)
  pngfile <- paste(directory, paste(toString(subid), '_wear_detection.png', sep=''), sep="/")
  plotWearGraph(acc, pngfile)
}

