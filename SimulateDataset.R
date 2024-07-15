### This script generates data that simulates n-days of accelerometer usage,
### where the weartime is from the GOALS study and the nonwear is from the Sleep study.
### Each simulated sample represents 1-day accelerometry readings for a simulated subject. 


### Main program

## Clear environment
rm(list = ls())

## Import libraries and related scripts
library(lubridate)
library(stringr)
# Load necessary libraries
library(dplyr)
library(tidyr)
source('WeartimeDetection.R')



## Global variables, including paths and simulation settings

# location of the goals data
goals_directory <- 'path/to/data_GOAL'
goals_accfilename_ending <- c('t00c11secDataTable.csv','t00c01secDataTable.csv','t00c21secDataTable.csv')

directory <-"path/to/output"

subjectids_goal = c(list_of_IDs)

files_goal <- unlist(lapply(subjectids_goal, function(subid) {
  pattern <- paste(toString(subid), "_weartime.csv", sep = "")
  list.files(path = directory, pattern = pattern, full.names = TRUE)
}))

# location of the sleep subject data
sleep_directory <- 'path/to/data_Sleep'
subjectids_sleep = c(list_of_IDs)
# List all files in the sleep directory and filter files that match the expected patterns
files_sleep <- list.files(path = sleep_directory, full.names = TRUE, pattern = "\\.csv$")

#apply choi and troiano to goal dataset
for (subid in subjectids_goal) {

  for (ending in goals_accfilename_ending) {
    # Generate the file path
    file_path <- paste(goals_directory, paste(toString(subid), ending, sep=''), sep = '/')
    
    # Check if the file exists
    if (file.exists(file_path)) {
      # If the file exists, assign the file path to goal_file and break out of the loop
      actigraphfile <- file_path
      break
    }
  }
  # generate weartime and actigraph file names based on subject id
  #actigraphfile <- paste(goals_directory, paste(toString(subid),goals_accfilename_ending, sep=''), sep = '/')
  print(actigraphfile)
  weartimefile <- paste(directory, paste(toString(subid), '_weartime.csv', sep=''), sep="/")

  # run weartime detection
  
  acc <-runWeartimeDetection(actigraphfile, directory,subid, usefixedchoi=TRUE, frame = 90, streamFrame = 30,numskiplines=10, window = 60, tol = 2, tol_upper=100, nci = FALSE, days_distinct = TRUE)
  acc$ID <- subid
  # write the weartime
  write.table(acc, file = weartimefile, sep = ",", row.names = F, col.name = T)
}


# Step 1: approach 1 Create a pool for each day 

daily_pool <- bind_rows(lapply(files_goal, function(file) {
  print(file)
  patient_data <- read.csv(file)
  # Extract ID from the file name
  id <- sub("_weartime.csv", "", basename(file))  # Adjust this based on your actual naming pattern
  patient_data %>%
    group_by(ID = id, Date = as.Date(Date, format = "%m/%d/%Y")) %>%
    filter(
      all(weartimeChoiVM == 1) & all(weartimeNhanesVM == 1) &&  # All timestamps must meet this condition
        n() == 86400  # Total count of timestamps for the day should be 86,400
    ) %>%
    ungroup()
})) %>%
  ungroup()
# Write daily pool as a CSV file
write.csv(daily_pool, file = file.path(directory, "daily_pool.csv"), row.names = FALSE)

# Step2.1 Define a function to split and assign to pools based on 15-minute intervals
split_and_create_15min <- function(patient_data) {
  patient_data %>%
    mutate(DateTime = as.POSIXct(paste(Date, Time), format = "%m/%d/%Y %H:%M:%S")) %>%
    mutate(interval = cut(DateTime, breaks = "15 min")) %>%
    group_by(id, Date, interval) %>%
    filter(
      all(weartimeTruth == 0)  # Condition for wear time
    ) %>%
    ungroup()
}

renamecolumn <- function(df, old_name, new_name) {
  if (old_name %in% names(df)) {
    names(df)[names(df) == old_name] <- new_name
    return(df)
  } else {
    warning(paste("Column", old_name, "does not exist in the dataframe. No changes made."))
    return(df)
  }
}

# Apply the function to each file and bind rows
nonwear_pool <- bind_rows(lapply(files_sleep, function(file) {
  patient_data <- read.csv(file)
  ts <- seq(c(ISOdate(2017,1,1)), by = "sec", length.out = nrow(patient_data))
  patient_data$Date <- strftime(ts, format="%m/%d/%Y")
  patient_data$Time <- strftime(ts, format="%H:%M:%S %p")
  
  # rename columns
  patient_data <- renamecolumn(patient_data, 'waistwaist_axis1', 'Axis1')
  patient_data <- renamecolumn(patient_data, 'waistwaist_axis2', 'Axis2')
  patient_data <- renamecolumn(patient_data, 'waistwaist_axis3', 'Axis3')
  patient_data <- renamecolumn(patient_data, 'waistwaist_steps', 'Steps')
  patient_data <- renamecolumn(patient_data, 'waistwaist_lux', 'Lux')
  patient_data <- renamecolumn(patient_data, 'waistwaist_incl_off', 'Inclinometer.Off')
  patient_data <- renamecolumn(patient_data, 'waistwaist_incl_stand', 'Inclinometer.Standing')
  patient_data <- renamecolumn(patient_data, 'waistwaist_incl_sit', 'Inclinometer.Sitting')
  patient_data <- renamecolumn(patient_data, 'waistwaist_incl_lying', 'Inclinometer.Lying')
  patient_data <- renamecolumn(patient_data, 'waistwaist_vm', 'Vector.Magnitude')
  split_and_create_15min(patient_data)
})) %>%
  ungroup()

# Write daily interval pool as a CSV file
write.csv(nonwear_pool, file = file.path(directory, "nonwear_pool_15_min.csv"), row.names = FALSE)

# Step2.2 Define a function to split and assign to pools based on 60-minute intervals
split_and_create_60min <- function(patient_data) {
  patient_data %>%
    mutate(DateTime = as.POSIXct(paste(Date, Time), format = "%m/%d/%Y %H:%M:%S")) %>%
    mutate(interval = cut(DateTime, breaks = "60 min")) %>%
    group_by(id, Date, interval) %>%
    filter(
      all(weartimeTruth == 0)  # Condition for wear time
    ) %>%
    ungroup()
}

# Apply the function to each file and bind rows
nonwear_pool <- bind_rows(lapply(files_sleep, function(file) {
  patient_data <- read.csv(file)
  ts <- seq(c(ISOdate(2017,1,1)), by = "sec", length.out = nrow(patient_data))
  patient_data$Date <- strftime(ts, format="%m/%d/%Y")
  patient_data$Time <- strftime(ts, format="%H:%M:%S %p")
  
  # rename columns
  patient_data <- renamecolumn(patient_data, 'waistwaist_axis1', 'Axis1')
  patient_data <- renamecolumn(patient_data, 'waistwaist_axis2', 'Axis2')
  patient_data <- renamecolumn(patient_data, 'waistwaist_axis3', 'Axis3')
  patient_data <- renamecolumn(patient_data, 'waistwaist_steps', 'Steps')
  patient_data <- renamecolumn(patient_data, 'waistwaist_lux', 'Lux')
  patient_data <- renamecolumn(patient_data, 'waistwaist_incl_off', 'Inclinometer.Off')
  patient_data <- renamecolumn(patient_data, 'waistwaist_incl_stand', 'Inclinometer.Standing')
  patient_data <- renamecolumn(patient_data, 'waistwaist_incl_sit', 'Inclinometer.Sitting')
  patient_data <- renamecolumn(patient_data, 'waistwaist_incl_lying', 'Inclinometer.Lying')
  patient_data <- renamecolumn(patient_data, 'waistwaist_vm', 'Vector.Magnitude')
  split_and_create_60min(patient_data)
})) %>%
  ungroup()

# Write daily interval pool as a CSV file
write.csv(nonwear_pool, file = file.path(directory, "nonwear_pool_1_hour.csv"), row.names = FALSE)


# Step 3: Create multiple simulations with random samples
#######################################################################

# Function to create a simulation dataset for daily pool option
create_daily_pool_simulation <- function(pool_data, n_days) {
  
  # Initialize an empty data frame to store simulation data
  simulation_data <- data.frame()
  #start_datetime <- as.POSIXct("2023-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%S")
  
  # Loop through each day
  for (day_index in 1:n_days) {
    
    # Randomly select one day
  
    sampled_rows <- pool_data %>%
      sample_n(size = 1, replace = TRUE)
    
    # Filter the entire pool data for each sampled day
    day <- pool_data %>%
      filter(ID %in% sampled_rows$ID, Date %in% sampled_rows$Date)
    day$DateTime <- as.POSIXct(paste(day$Date, day$Time), format = "%m/%d/%Y %H:%M:%S")
    # Append the interval data to the simulation dataset
    simulation_data <- bind_rows(simulation_data, day)
  }
  #simulation_data$DateTime_new <- start_datetime + (0:(nrow(simulation_data)-1))  # Increase by 1 sec per row
  return(simulation_data)
}

generate_nonwear_period_less_than_one_hour <- function(nonwear_pool) {
  # Initialize an empty data frame to store the nonwear period
  nonwear_period <- data.frame()
  
  # Define parameters
  threshold <- 1  # Probability threshold
  max_duration <- 60 * 60  # Maximum duration in seconds (1 hour)
  interval_length <- 15 * 60  # Length of each 15-minute interval in seconds
  
  # Continue adding 15-minute intervals until the length of the period exceeds 1 hour
  while (nrow(nonwear_period) < max_duration) {
    # Generate a random number between 0 and 1 for probability
    probability <- runif(1)
    
    # If the probability is less than the threshold, add a 15-minute interval to the nonwear period
    if (probability < threshold) {
      # Sample a random interval from nonwear_pool with replacement 
      nonwear_interval <- nonwear_pool %>%
        distinct(id, interval) %>%
        sample_n(size = 1, replace = TRUE)
      
      # Filter nonwear_pool to get the data for the sampled interval
      nonwear_data <- nonwear_pool %>%
        filter(id %in% nonwear_interval$id, interval %in% nonwear_interval$interval)
      nonwear_data$weartimeTruth <- 0
      # Add the sampled interval to the nonwear period
      nonwear_period <- rbind(nonwear_period, nonwear_data)
      threshold <- threshold * 0.8
      
    } else {
      # Break the loop if the probability is greater than or equal to the threshold
      break
    }
  }
  
  return(nonwear_period)
}

generate_nonwear_period_greater_than_one_hour <- function(nonwear_pool) {
  # Initialize an empty data frame to store the nonwear period
  nonwear_period <- data.frame()
  
  # Define parameters
  threshold <- 1  # Probability threshold
  max_duration <- 24 * 60 * 60  # Maximum duration in seconds (24 hours)
  interval_length <- 60 * 60  # Length of each 1-hour interval in seconds
  
  # Continue adding 1-hour intervals until the length of the period exceeds 24 hours
  while (nrow(nonwear_period) < max_duration) {
    # Generate a random number between 0 and 1 for probability
    probability <- runif(1)
    
    # If the probability is less than the threshold, add a 1hour interval to the nonwear period
    if (probability < threshold) {
      # Sample a random interval from nonwear_pool with replacement 
      nonwear_interval <- nonwear_pool %>%
        distinct(id, interval) %>%
        sample_n(size = 1, replace = TRUE)
      
      # Filter nonwear_pool to get the data for the sampled interval
      nonwear_data <- nonwear_pool %>%
        filter(id %in% nonwear_interval$id, interval %in% nonwear_interval$interval)
      nonwear_data$weartimeTruth <- 0
      # Add the sampled interval to the nonwear period
      nonwear_period <- rbind(nonwear_period, nonwear_data)
      threshold <- threshold * 0.8
      
    } else {
      # Break the loop if the probability is greater than or equal to the threshold
      break
    }
  }
  return(nonwear_period)
}

replace_wear_with_nonwear <- function(day_data_final,day_data, nonwear_periods, day_num) {
  for (j in 1:2) {
    nonwear_data <- nonwear_periods[[j]]
    nonwear_data <- subset(nonwear_data, select = c("id", "DateTime", "Axis1", "Axis2", "Axis3","Vector.Magnitude","weartimeTruth"))
    names(nonwear_data) <- c("ID", "DateTime", "Axis1", "Axis2", "Axis3","Vector.Magnitude","weartimeTruth")  
    
    if (j == 1) {
      day_data[1:nrow(nonwear_data), ] <- nonwear_data
    } else {
      day_data[(nrow(day_data)-nrow(nonwear_data) + 1):nrow(day_data), ] <- nonwear_data
    }
  }
  if (day_num == 1) {
    day_data_final <- day_data
  } else {
    day_data_final <- rbind(day_data_final, day_data)
  }
  return(day_data_final)
}


# Main function to simulate dataset
simulate_dataset <- function(num_simulations, daily_pool, n_days, nonwear_pool_15_min, nonwear_pool_1_hour) {
  simulation_list <- list()
  
  for (sim_num in 1:num_simulations) {
    print(sim_num)
    # Set seed for reproducibility
    set.seed(sim_num)
    # Create a simulation dataset for daily_pool option
    daily_pool_simulation <- create_daily_pool_simulation(daily_pool, n_days)
    print('daily_pool_simulation is ready')
    # Generate nonwear periods
    nonwear_periods <- list()
    for (i in 1:2) { # Generate 2 nonwear periods
      threshold <- runif(1)
      # If the threshold is less than 0.10, use generate_nonwear_period_less_than_one_hour, otherwise use generate_nonwear_period_greater_than_one_hour 
      if (threshold < 0.10) {
        nonwear_period <- generate_nonwear_period_less_than_one_hour(nonwear_pool_15_min)
      } else {
        nonwear_period <- generate_nonwear_period_greater_than_one_hour(nonwear_pool_1_hour)
      }
      nonwear_periods[[i]] <- nonwear_period
    }
    print('2 nonwear durations are ready')
    day_data_final <- data.frame()
    probability1 <- runif(1)
    # If the probability1 is less than 0.20, add 24 hours knocked out days with some probability 
    if (probability1 >= 0.20) {
      for (day_num in 1:n_days) {
        day_data <- daily_pool_simulation[((day_num - 1) * 86400 + 1):(day_num * 86400), ]
        print(nrow(day_data))
        day_data <- subset(day_data, select = c("ID", "DateTime", "Axis1", "Axis2", "Axis3", "Vector.Magnitude"))  # Adjust column names as needed
        day_data$weartimeTruth <- 1
        day_data_final <- replace_wear_with_nonwear(day_data_final,day_data, nonwear_periods, day_num)
      }
    } else {
      for (day_num in 1:n_days) {
        probability2 <- runif(1)
        if (probability2 < 0.40) {
          nonwear_interval <- nonwear_pool_1_hour %>%
            group_by(id, interval) %>%
            filter(n() == 3600) %>%
            distinct(id, interval) %>%
            ungroup() %>%
            sample_n(size = 24, replace = FALSE)
          day_data <- nonwear_pool_1_hour %>%
            semi_join(nonwear_interval, by = c("id", "interval"))
          day_data$weartimeTruth <- 0
          day_data <- subset(day_data, select = c("id", "DateTime", "Axis1", "Axis2", "Axis3", "Vector.Magnitude", "weartimeTruth"))  # Adjust column names as needed
          names(day_data) <- c("ID", "DateTime", "Axis1", "Axis2", "Axis3", "Vector.Magnitude", "weartimeTruth")  
          print(nrow(day_data))
          if (day_num == 1) {
            day_data_final <- day_data
          } else {
            day_data_final <- rbind(day_data_final, day_data)
          }
        } else {
          day_data <- daily_pool_simulation[((day_num - 1) * 86400 + 1):(day_num * 86400), ]
          day_data <- subset(day_data, select = c("ID", "DateTime", "Axis1", "Axis2", "Axis3", "Vector.Magnitude"))
          day_data$weartimeTruth <- 1
          day_data_final <- replace_wear_with_nonwear(day_data_final,day_data, nonwear_periods, day_num)
        }
      }
    }
    start_datetime <- as.POSIXct("2023-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%S")
    day_data_final$DateTime_new <- start_datetime + (0:(nrow(day_data_final)-1))  # Increase by 1 sec per row
    day_data_final$ID_new <- sim_num
    # Save the data as CSV
    simulation_name <- paste0("daily_pool_simulation_", sim_num, ".csv")
    simulationfile <- file.path(simulation_directory, simulation_name)
    write.csv(day_data_final, file = simulationfile, row.names = FALSE)
    simulation_list[[sim_num]] <- day_data_final
  }
  
  return(simulation_list)
}


# folder for the simulated data
simulation_directory <- 'path/to/data_simulation'

# Define the number of simulations and number of days
num_simulations <- 50  # Adjust as needed
n_days <- 1  # Number of days to simulate

sim_data<-simulate_dataset(num_simulations, daily_pool, n_days, nonwear_pool_15_min, nonwear_pool_1_hour) 


