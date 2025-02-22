### This script generates data that simulates 24-hours of accelerometer usage,
### where the weartime is from the GOALS and Sleep study and the nonwear is from the Sleep study.
### Each simulated sample represents 24-hours accelerometry readings for a simulated subject. 

# Clear environment
rm(list = ls())

# Import libraries
library(dplyr)
library(readr)

# location of the goals data
directory <- 'path/to/data_GOAL'
output_directory <-"path/to/output"

subjectids_goal = c(list_of_IDs)

# List files based on subject IDs
files_goal <- unlist(lapply(subjectids_goal, function(subid) {
  pattern <- paste(toString(subid), "_weartime.csv", sep = "")
  list.files(path = directory, pattern = pattern, full.names = TRUE)
}))

# location of the sleep subject data
sleep_directory <- 'path/to/data_Sleep'
subjectids_sleep = c(list_of_IDs)
# List all files in the sleep directory and filter files that match the expected patterns
files_sleep <- list.files(path = sleep_directory, full.names = TRUE, pattern = "\\.csv$")

################################################################################
#Create non-sedentary pool

process_data <- function(file) {
  print(paste("Processing file:", file))
  
  patient_data <- read_csv(file, show_col_types = FALSE)
  print(paste("Total rows loaded:", nrow(patient_data)))
  patient_data$TimeStamp <- as.POSIXct(patient_data$TimeStamp, format = "%Y-%m-%d %H:%M:%S")
  
  patient_data <- patient_data %>%
    mutate(row_id = row_number(), 
           interval = floor(as.numeric(difftime(TimeStamp, min(TimeStamp), units = "secs")) / 15))
  
  print(paste("Intervals created:", length(unique(patient_data$interval))))
  
  patient_data_processed <- patient_data %>%
    group_by(interval) %>%
    summarize(
      AxisVM_sum = sum(Vector.Magnitude, na.rm = TRUE), 
      TimeStamp = first(TimeStamp),
      StartRow = first(row_id),
      EndRow = last(row_id),
      .groups = 'drop'
    ) %>%
    mutate(group = cumsum(AxisVM_sum <= 180)) %>%
    group_by(group) %>%
    filter(AxisVM_sum > 180) %>%
    mutate(duration = n() * 15) %>%
    filter(duration <= 14400 & duration >= 900) %>%
    summarize(
      StartTime = first(TimeStamp),
      EndTime = last(TimeStamp),
      TotalDuration = sum(duration),
      StartRow = first(StartRow),
      EndRow = last(EndRow),
      .groups = 'drop'
    )
  
  print(paste("Segments processed:", nrow(patient_data_processed)))
  
  segments_list <- list()
  for (i in seq_len(nrow(patient_data_processed))) {
    segment <- patient_data %>%
      filter(row_id >= patient_data_processed$StartRow[i] & 
               row_id <= patient_data_processed$EndRow[i])
    
    if (nrow(segment) > 0) {
      segments_list[[length(segments_list) + 1]] <- segment
    } else {
      print(paste("No data found for rows:", patient_data_processed$StartRow[i], "to", patient_data_processed$EndRow[i]))
    }
  }
  
  return(segments_list)
}

# Apply function to each file and store results in a list of lists
nonsedentary_pools <- lapply(files_goal, process_data)

saveRDS(nonsedentary_pools, file.path(output_directory, "nonsedentary_pools_vm.rds"))

################################################################################
#Create sedentary and nonwear pools

process_and_extract_segments <- function(file) {
  patient_data <- read.csv(file)
  patient_data$DateTime <- as.POSIXct(paste(patient_data$Date, patient_data$Time), format = "%m/%d/%Y %H:%M:%S")
  
  # Generate non-overlapping intervals
  current_time <- min(patient_data$DateTime)
  intervals <- data.frame(StartDateTime = numeric(), EndDateTime = numeric())
  
  while(current_time < max(patient_data$DateTime)) {
    interval_length <- sample(900:14400, 1)  # Duration between 15 min and 4 hours
    end_time <- min(current_time + interval_length, max(patient_data$DateTime))
    
    intervals <- rbind(intervals, data.frame(StartDateTime = current_time, EndDateTime = end_time))
    current_time <- end_time + 1  # Ensure non-overlapping intervals
  }
  
  # Assign data points to intervals
  patient_data$IntervalID <- findInterval(patient_data$DateTime, intervals$EndDateTime, all.inside = TRUE)
  
  # Extract actual data for each interval
  segments_list <- list()
  for (i in 1:nrow(intervals)) {
    segment <- patient_data %>%
      filter(DateTime >= intervals$StartDateTime[i] & DateTime <= intervals$EndDateTime[i])
    
    # Determine wear time truth
    if (all(segment$weartimeTruth == 0)) {
      wear_time_label <- "nonwear"
    } else if (all(segment$weartimeTruth == 1)) {
      wear_time_label <- "sedentary"
    } else {
      wear_time_label <- "mixed"
    }
    
    # Store segment if it meets criteria
    if (wear_time_label %in% c("nonwear", "sedentary")) {
      segments_list[[length(segments_list) + 1]] <- list(
        Data = segment,
        WearTimeLabel = wear_time_label
      )
    }
  }
  
  return(segments_list)
}

# Apply the function to each file
all_segments <- lapply(files_sleep, process_and_extract_segments)

saveRDS(all_segments, file.path(output_directory, "nonwear_sedentary_pool.rds"))

################################################################################
# Extract long/short sedentary, non-sedentary and nonwear segments 

# Function to extract relevant segments from a given pool
extract_segments <- function(pool, label, min_non_zero_count=NULL,max_non_zero_count=NULL,min_duration=NULL, max_duration=NULL) {
  segments <- list()
  
  for (i in seq_along(pool)) {
    for (j in seq_along(pool[[i]])) {
      segment <- nonwear_sedentary_pools[[i]][[j]]
      # Check the label and optional duration conditions
      if (segment$WearTimeLabel == label &&
          (sum(segment$Data["Vector.Magnitude"] != 0) / nrow(segment$Data["Vector.Magnitude"]))* 100<max_non_zero_count &&
          (sum(segment$Data["Vector.Magnitude"] != 0) / nrow(segment$Data["Vector.Magnitude"]))* 100>=min_non_zero_count &&
          (nrow(segment$Data) >= min_duration) &&
          (nrow(segment$Data)<= max_duration)) {
        segments[[length(segments) + 1]] <- segment$Data
      }
    }
  }
  
  return(segments)
}

extract_segments_nonwear <- function(pool, label,min_duration=NULL, max_duration=NULL) {
  segments <- list()
  
  for (i in seq_along(pool)) {
    for (j in seq_along(pool[[i]])) {
      segment <- nonwear_sedentary_pools[[i]][[j]]
      # Check the label and optional duration conditions
      if (segment$WearTimeLabel == label &&
          (nrow(segment$Data) >= min_duration) &&
          (nrow(segment$Data)<= max_duration)) {
        segments[[length(segments) + 1]] <- segment$Data
      }
    }
  }
  
  return(segments)
}

extract_nonsedentary_segments <- function(pool,min_duration=NULL, max_duration=NULL) {
  segments <- list()
  for (i in seq_along(pool)) {
    for (j in seq_along(pool[[i]])) {
      segment <- pool[[i]][[j]]
      # Check the label and optional duration conditions
      if ((nrow(segment) >= min_duration) &&
          (nrow(segment)<= max_duration)) {
        segments[[length(segments) + 1]] <- segment
      }
    }
  }
  return(segments)
}

# Extract segments based on labels and conditions
sedentary_segments_long <- extract_segments(nonwear_sedentary_pools , "sedentary", min_non_zero_count=2,max_non_zero_count=100.1, min_duration = 3600, max_duration = 14400)
sedentary_segments_short <- extract_segments(nonwear_sedentary_pools , "sedentary", min_non_zero_count=2,max_non_zero_count=100.1, min_duration = 900, max_duration = 3600)
sleep_segments_long <- extract_segments(nonwear_sedentary_pools , "sedentary", min_non_zero_count=0,max_non_zero_count=2, min_duration = 3600, max_duration = 14400)
sleep_segments_short <- extract_segments(nonwear_sedentary_pools , "sedentary", min_non_zero_count=0,max_non_zero_count=2, min_duration = 900, max_duration = 3600)
nonsedentary_segments <- extract_nonsedentary_segments(nonsedentary_pools,min_duration = 900, max_duration = 14400)
nonwear_segments_short <- extract_segments_nonwear(nonwear_sedentary_pools , "nonwear", min_duration = 900, max_duration = 3600)  # Short nonwear
nonwear_segments_long <- extract_segments_nonwear(nonwear_sedentary_pools , "nonwear", min_duration = 3600,  max_duration = 14400)  # Long nonwear

################################################################################
# Simulate data

simulate_24hr_data <- function(seed) {
  set.seed(seed)
  total_duration <- 0  # Total accumulated time in seconds
  combined_data <- list()  # To store the combined simulated data
  
  # Step 1: Simulate daytime up to 16-18 hours (57600 seconds)
  while (total_duration < 64800) {
    random_step <- runif(1)
    
    if (random_step < 0.10) {
      # Add a non-sedentary segment
      nonsedentary_sample_index <- sample(seq_along(nonsedentary_segments), 1)
      nonsedentary_sample <- nonsedentary_segments[[nonsedentary_sample_index]]
      nonsedentary_sample <- subset(nonsedentary_sample, select = c("ID", "TimeStamp", "Axis1", "Axis2", "Axis3", "Vector.Magnitude"))
      names(nonsedentary_sample) <- c("ID", "DateTime", "Axis1", "Axis2", "Axis3", "Vector.Magnitude", "weartimeTruth")  
      nonsedentary_sample$weartimeTruth <- 1 
      combined_data <- bind_rows(combined_data, nonsedentary_sample)
      total_duration <- total_duration + nrow(nonsedentary_sample)
      
    } else if (random_step >= 0.10 && random_step < 0.30) {
      # Add a sedentary short segment
      sedentary_short_sample_index <- sample(seq_along(sedentary_segments_short), 1)
      sedentary_short_sample <- sedentary_segments_short[[sedentary_short_sample_index]]
      sedentary_short_sample <- subset(sedentary_short_sample, select = c("id", "DateTime", "Axis1", "Axis2", "Axis3", "Vector.Magnitude", "weartimeTruth"))
      names(sedentary_short_sample) <- c("ID", "DateTime", "Axis1", "Axis2", "Axis3", "Vector.Magnitude", "weartimeTruth")  
      combined_data <- bind_rows(combined_data, sedentary_short_sample)
      total_duration <- total_duration + nrow(sedentary_short_sample)
      
    } else if (random_step >= 0.30 && random_step < 0.50) {
      # Add a sedentary long segment
      sedentary_long_sample_index <- sample(seq_along(sedentary_segments_long), 1)
      sedentary_long_sample <- sedentary_segments_long[[sedentary_long_sample_index]]
      sedentary_long_sample <- subset(sedentary_long_sample, select = c("id", "DateTime", "Axis1", "Axis2", "Axis3", "Vector.Magnitude", "weartimeTruth"))
      names(sedentary_long_sample) <- c("ID", "DateTime", "Axis1", "Axis2", "Axis3", "Vector.Magnitude", "weartimeTruth")  
      combined_data <- bind_rows(combined_data, sedentary_long_sample)
      total_duration <- total_duration + nrow(sedentary_long_sample)
      
    }else if (random_step >= 0.50 && random_step < 0.60) {
      # Add a nonwear short segment
      nonwear_short_sample_index <- sample(seq_along(nonwear_segments_short), 1)
      nonwear_short_sample <- nonwear_segments_short[[nonwear_short_sample_index]]
      nonwear_short_sample <- subset(nonwear_short_sample, select = c("id", "DateTime", "Axis1", "Axis2", "Axis3", "Vector.Magnitude", "weartimeTruth"))
      names(nonwear_short_sample) <- c("ID", "DateTime", "Axis1", "Axis2", "Axis3", "Vector.Magnitude", "weartimeTruth") 
      combined_data <- bind_rows(combined_data, nonwear_short_sample)
      total_duration <- total_duration + nrow(nonwear_short_sample)
      
    }else if (random_step >= 0.60 && random_step < 0.80) {
      # Add a nonwear long segment
      nonwear_long_sample_index <- sample(seq_along(nonwear_segments_long), 1)
      nonwear_long_sample <- nonwear_segments_long[[nonwear_long_sample_index]]
      nonwear_long_sample <- subset(nonwear_long_sample, select = c("id", "DateTime", "Axis1", "Axis2", "Axis3", "Vector.Magnitude", "weartimeTruth"))
      names(nonwear_long_sample) <- c("ID", "DateTime", "Axis1", "Axis2", "Axis3", "Vector.Magnitude", "weartimeTruth") 
      combined_data <- bind_rows(combined_data, nonwear_long_sample)
      total_duration <- total_duration + nrow(nonwear_long_sample)
      
    } else if (random_step >= 0.80 && random_step < 0.90) {
      # Add a short sleep segment (less common during the day)
      sleep_short_sample_index <- sample(seq_along(sleep_segments_short), 1)
      sleep_short_sample <- sleep_segments_short[[sleep_short_sample_index]]
      sleep_short_sample <- subset(sleep_short_sample, select = c("id", "DateTime", "Axis1", "Axis2", "Axis3", "Vector.Magnitude", "weartimeTruth"))
      names(sleep_short_sample) <- c("ID", "DateTime", "Axis1", "Axis2", "Axis3", "Vector.Magnitude", "weartimeTruth") 
      combined_data <- bind_rows(combined_data, sleep_short_sample)
      total_duration <- total_duration + nrow(sleep_short_sample)
    } else {
      # Add a long sleep segment (less common during the day)
      sleep_long_sample_index <- sample(seq_along(sleep_segments_long), 1)
      sleep_long_sample <- sleep_segments_long[[sleep_long_sample_index]]
      sleep_long_sample <- subset(sleep_long_sample, select = c("id", "DateTime", "Axis1", "Axis2", "Axis3", "Vector.Magnitude", "weartimeTruth"))
      names(sleep_long_sample) <- c("ID", "DateTime", "Axis1", "Axis2", "Axis3", "Vector.Magnitude", "weartimeTruth") 
      combined_data <- bind_rows(combined_data, sleep_long_sample)
      total_duration <- total_duration + nrow(sleep_long_sample)
    }
  }
  
  # Step 2: Simulate nighttime (the remaining time up to 24 hours)
  while (total_duration < 86400) {  # 24 hours
    random_step <- runif(1)
    
    if (random_step < 0.04) {
      # Add a non-sedentary segment (rare at night)
      nonsedentary_sample_index <- sample(seq_along(nonsedentary_segments), 1)
      nonsedentary_sample <- nonsedentary_segments[[nonsedentary_sample_index]]
      nonsedentary_sample <- subset(nonsedentary_sample, select = c("ID", "TimeStamp", "Axis1", "Axis2", "Axis3", "Vector.Magnitude"))
      names(nonsedentary_sample) <- c("ID", "DateTime", "Axis1", "Axis2", "Axis3", "Vector.Magnitude", "weartimeTruth")  
      nonsedentary_sample$weartimeTruth <- 1 
      combined_data <- bind_rows(combined_data, nonsedentary_sample)
      total_duration <- total_duration + nrow(nonsedentary_sample)
      
    } else if (random_step >= 0.04 && random_step < 0.20) {
      # Add a sedentary short segment
      sedentary_short_sample_index <- sample(seq_along(sedentary_segments_short), 1)
      sedentary_short_sample <- sedentary_segments_short[[sedentary_short_sample_index]]
      sedentary_short_sample <- subset(sedentary_short_sample, select = c("id", "DateTime", "Axis1", "Axis2", "Axis3", "Vector.Magnitude", "weartimeTruth"))
      names(sedentary_short_sample) <- c("ID", "DateTime", "Axis1", "Axis2", "Axis3", "Vector.Magnitude", "weartimeTruth")  
      combined_data <- bind_rows(combined_data, sedentary_short_sample)
      total_duration <- total_duration + nrow(sedentary_short_sample)
      
    } else if (random_step >= 0.20 && random_step < 0.36) {
      # Add a sedentary long segment
      sedentary_long_sample_index <- sample(seq_along(sedentary_segments_long), 1)
      sedentary_long_sample <- sedentary_segments_long[[sedentary_long_sample_index]]
      sedentary_long_sample <- subset(sedentary_long_sample, select = c("id", "DateTime", "Axis1", "Axis2", "Axis3", "Vector.Magnitude", "weartimeTruth"))
      names(sedentary_long_sample) <- c("ID", "DateTime", "Axis1", "Axis2", "Axis3", "Vector.Magnitude", "weartimeTruth")  
      combined_data <- bind_rows(combined_data, sedentary_long_sample)
      total_duration <- total_duration + nrow(sedentary_long_sample)
      
    }else if (random_step >= 0.36 && random_step < 0.46) {
      # Add a nonwear short segment
      nonwear_short_sample_index <- sample(seq_along(nonwear_segments_short), 1)
      nonwear_short_sample <- nonwear_segments_short[[nonwear_short_sample_index]]
      nonwear_short_sample <- subset(nonwear_short_sample, select = c("id", "DateTime", "Axis1", "Axis2", "Axis3", "Vector.Magnitude", "weartimeTruth"))
      names(nonwear_short_sample) <- c("ID", "DateTime", "Axis1", "Axis2", "Axis3", "Vector.Magnitude", "weartimeTruth") 
      combined_data <- bind_rows(combined_data, nonwear_short_sample)
      total_duration <- total_duration + nrow(nonwear_short_sample)
      
    }else if (random_step >= 0.46 && random_step < 0.66) {
      # Add a nonwear long segment
      nonwear_long_sample_index <- sample(seq_along(nonwear_segments_long), 1)
      nonwear_long_sample <- nonwear_segments_long[[nonwear_long_sample_index]]
      nonwear_long_sample <- subset(nonwear_long_sample, select = c("id", "DateTime", "Axis1", "Axis2", "Axis3", "Vector.Magnitude", "weartimeTruth"))
      names(nonwear_long_sample) <- c("ID", "DateTime", "Axis1", "Axis2", "Axis3", "Vector.Magnitude", "weartimeTruth") 
      combined_data <- bind_rows(combined_data, nonwear_long_sample)
      total_duration <- total_duration + nrow(nonwear_long_sample)
      
    }else if (random_step >= 0.66 && random_step < 0.83) {
      # Add a short sleep segment (less common during the day)
      sleep_short_sample_index <- sample(seq_along(sleep_segments_short), 1)
      sleep_short_sample <- sleep_segments_short[[sleep_short_sample_index]]
      sleep_short_sample <- subset(sleep_short_sample, select = c("id", "DateTime", "Axis1", "Axis2", "Axis3", "Vector.Magnitude", "weartimeTruth"))
      names(sleep_short_sample) <- c("ID", "DateTime", "Axis1", "Axis2", "Axis3", "Vector.Magnitude", "weartimeTruth") 
      combined_data <- bind_rows(combined_data, sleep_short_sample)
      total_duration <- total_duration + nrow(sleep_short_sample)
      
    } else {
      # Add a long sleep segment (less common during the day)
      sleep_long_sample_index <- sample(seq_along(sleep_segments_long), 1)
      sleep_long_sample <- sleep_segments_long[[sleep_long_sample_index]]
      sleep_long_sample <- subset(sleep_long_sample, select = c("id", "DateTime", "Axis1", "Axis2", "Axis3", "Vector.Magnitude", "weartimeTruth"))
      names(sleep_long_sample) <- c("ID", "DateTime", "Axis1", "Axis2", "Axis3", "Vector.Magnitude", "weartimeTruth") 
      combined_data <- bind_rows(combined_data, sleep_long_sample)
      total_duration <- total_duration + nrow(sleep_long_sample)
    }
    
    # If the total duration exceeds 24 hours, trim the last segment
    if (total_duration > 86400) {
      combined_data <- combined_data[1:86400, ]
      total_duration <- 86400  # Set to exactly 24 hours
    }
  }
  
  # Add a DateTime column to the combined data for tracking time progression
  combined_data$DateTime <- seq(from = ymd_hms("2023-01-01 00:00:00"), by = "secs", length.out = nrow(combined_data))
  
  return(combined_data)
}

# Function to simulate multiple datasets and save them
simulate_and_save_datasets <- function(n, output_directory) {
  for (i in 1:n) {
    simulated_data <- simulate_24hr_data(seed=i)
    file_name <- sprintf("simulated_data_%02d.csv", i)
    write.csv(simulated_data, file.path(output_directory, file_name), row.names = FALSE)
  }
}

n_datasets <- 100

simulation_directory <-'path/to/simulation_data'
# Generate and save datasets
simulate_and_save_datasets(n_datasets, simulation_directory)
