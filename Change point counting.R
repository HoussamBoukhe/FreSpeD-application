library(tibble)
library(dplyr)
library(readr)

setwd("C:/Users/labo-maths/Desktop/FreSpeD-application")

#sum_subject <- read_rds("Output/windowLen125-Welch1.rds")

get_cp_count <- function(rds_file_path){
  
  sum_subject <- read_rds(rds_file_path)
  freq_band_names <- rownames(sum_subject$cp_sum[[1]]$tv_estimates$AF7)
  num_fb <- length(freq_band_names)
  
  # Create empty data frame
  cp_count_sum <- data.frame(subject = character(),
                             region = integer(),
                             epoch = integer(),
                             eyeState = integer()
  )
  for (freq in freq_band_names) {
    cp_count_sum[[paste0("nCp_", freq)]] <- integer()
  }
  
  n <- dim(sum_subject)[1]
  for (i in 1:n){ 
    tmp_data <- sum_subject[i,]
    tmp_subject <- tmp_data$subject
    tmp_region <- tmp_data$region
    tmp_epoch <- tmp_data$epoch
    tmp_eyeState <- tmp_data$eyeState
    tmp_cusum_values <- tmp_data$cp_sum[[1]]$cusum_values
    tmp_count <- rep(0, num_fb)
    
    for (j in 1:length(tmp_cusum_values)){
      tmp_count <- tmp_count + rowSums(tmp_cusum_values[[j]]>0)
    }
    
    if ((i %% 1000) == 0){print(i) } # ...
    
    # Create a new row as a named list or data frame
    new_row <- data.frame(
      subject = tmp_subject,
      region = tmp_region,
      epoch = tmp_epoch,
      eyeState = tmp_eyeState
    )
    
    # Add the frequency columns
    for (k in 1:length(freq_band_names)) {
      new_row[[paste0("nCp_", freq_band_names[k])]] <- tmp_count[k]
    }
    
    # Bind the properly structured row
    cp_count_sum <- rbind(cp_count_sum, new_row)
  }
  
  return(cp_count_sum)
}

sum1 <- get_cp_count("Output/windowLen125-Welch1.rds")
sum2 <- get_cp_count("Output/windowLen250-Welch2.rds")
sum3 <- get_cp_count("Output/windowLen500-Welch4.rds")

write.csv(sum1,"Output/windowLen125-Welch1.csv")
write.csv(sum2,"Output/windowLen250-Welch2.csv")
write.csv(sum3,"Output/windowLen500-Welch4.csv")
