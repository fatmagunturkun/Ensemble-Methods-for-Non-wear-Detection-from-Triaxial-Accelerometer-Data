# get confusion matrix of different predictors
# wt is the data frame of weartime predictions
# gt is the ground truth column

getConfusionMatrix <- function(wt, gt)
{
    cols <- colnames(wt)
    m <- ncol(wt)
    n <- nrow(wt)
    confusion <-  data.frame(subid = 0)
    for (i in (1:m)) {
        for (v in unique(gt)) {
            retrieved <- sum(wt[,i]==v)
            if (retrieved > 0) {
                confusion[, paste(cols[i],"_precision",toString(v), sep ="")] <- sum((wt[,i]==v) & (gt==v))/retrieved
            } else {
                confusion[, paste(cols[i],"_precision",toString(v), sep ="")] <- 1
            }
            relevant <- sum(gt==v)
            if (relevant > 0) {
                confusion[, paste(cols[i],"_recall",toString(v), sep ="")] <- sum((wt[,i]==v) & (gt==v))/relevant
            } else {
                confusion[, paste(cols[i],"_recall",toString(v), sep ="")] <- 1
            }
            negative <- sum(gt!=v)
            if (negative > 0) {
                confusion[, paste(cols[i],"_specificity",toString(v), sep ="")] <- sum((wt[,i]!=v) & (gt!=v))/negative
            } else {
                confusion[, paste(cols[i],"_specificity",toString(v), sep ="")] <- 1
            }
        }
    }
    return (confusion)
}

### Main program

library(lubridate)
library(stringr)
