library(tibble)
library(dplyr)
library(readr)

setwd("C:/Users/labo-maths/Desktop/FreSpeD-application")

#sum_subject <- read_rds("Output/windowLen125-Welch1.rds")

library(tibble)
library(dplyr)
library(readr)

get_change_point_loc <- function(rds_file_path){
  sum_subject <- read_rds(rds_file_path)
  freq_band_names <- rownames(sum_subject$cp_sum[[1]]$tv_estimates$AF7)
  num_fb <- length(freq_band_names)
  
  # Create empty data frame
  cp_count_sum <- data.frame(subject = character(),
                             region = integer(),
                             epoch = integer(),
                             eyeState = integer())
  
  # Initialize columns for each frequency band
  for (freq in freq_band_names) {
    cp_count_sum[[paste0("CpLoc_", freq)]] <- list()
  }
  
  n <- dim(sum_subject)[1]
  for (i in 1:n){
    tmp_data <- sum_subject[i,]
    tmp_subject <- tmp_data$subject 
    tmp_region <- tmp_data$region
    tmp_epoch <- tmp_data$epoch
    tmp_eyeState <- tmp_data$eyeState
    
    # Create a named list for the frequency bands
    tmp_loc_all_freq <- list()
    
    tmp_cusum_values <- tmp_data$cp_sum[[1]]$cusum_values    
    for (fb in 1:num_fb){
      freq_name <- freq_band_names[fb]
      tmp_fb_loc <- c()
      for (j in 1:length(tmp_cusum_values)){
        tmp_loc <- which(tmp_cusum_values[[j]][fb,] > 0)
        if (length(tmp_loc) > 0){
          tmp_fb_loc <- c(tmp_fb_loc, tmp_loc)
        }
      }
      tmp_fb_loc <- sort(tmp_fb_loc)
      # Store with the correct frequency name
      tmp_loc_all_freq[[paste0("CpLoc_", freq_name)]] <- list(tmp_fb_loc)
    }
    
    # Add a new row to the dataframe
    new_row <- data.frame(
      subject = tmp_subject,
      region = tmp_region,
      epoch = tmp_epoch,
      eyeState = tmp_eyeState
    )
    
    # Add frequency band change point locations
    for (freq in freq_band_names) {
      new_row[[paste0("CpLoc_", freq)]] <- tmp_loc_all_freq[[paste0("CpLoc_", freq)]]
    }
    
    # Append the new row
    cp_count_sum <- bind_rows(cp_count_sum, new_row)
  }
  
  return(cp_count_sum)
}



# Function to convert list columns to JSON strings for CSV export
prepare_dataframe_for_csv <- function(df) {
  # Load jsonlite library
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    install.packages("jsonlite")
  }
  library(jsonlite)
  
  # Create a copy of the dataframe
  df_for_csv <- df
  
  # Identify list columns
  list_cols <- sapply(df, is.list)
  
  # Convert list columns to JSON strings
  for (col in names(df)[list_cols]) {
    df_for_csv[[col]] <- sapply(df[[col]], function(x) {
      if (length(x) == 0) {
        return("[]")
      } else {
        return(toJSON(x))
      }
    })
  }
  
  return(df_for_csv)
}

windowLen <- "500"
M_welch <- "4"

rds_file_path <- paste("Output/windowLen",windowLen,"-Welch"
                       ,M_welch,".rds",sep = "")
cp_loc <- get_change_point_loc(rds_file_path)
tmp_for_csv <- prepare_dataframe_for_csv(cp_loc)
write.csv(tmp_for_csv, paste("Output/location-windowLen",windowLen,"-Welch",
                             M_welch,".csv",sep=""), row.names = FALSE)
