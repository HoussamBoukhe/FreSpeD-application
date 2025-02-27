#'@description
#'In this file, the FreSpeD algorithm is applied to each EEG recording 
#'for the 204 subject with open and closed eyes.
#'
#' @output The output is a *.rds file that contains a dataframe with 
#' the following columns: 
#' 
#' @column1 subject, charachter. Contains a subject unique id
#' @column2 region, integer in c(1:7). The algorithm is applied only 
#' to this electrodes
#' @column3 epoch, integer in c(1:8). Each recording has 8 minutes, 
#' the algorithm is applied to an interval of 1 minute.
#' @column4 eyeState, integer in c(0,1). 0 for eye closed and 1 for open
#' @column5 cp_sum, list. A list containing the output of FreSpeD  
#' applied to the EEG data, with the specified four columns.



setwd("C:/Users/labo-maths/Desktop/FreSpeD-application") # CHANGE with the location of your files

source("FreSpeD.r")
# library(remotes)
# remotes::install_github("craddm/eegUtils")
library(eegUtils) # installed from github repository to read .set files 
library(tibble)
library(dplyr)
library(readr) #save the output of the fresped algorithm 

################################ Regions of interest Roi_*
Roi_1 <- c("AF7","F7","F5","F3","FT7","FC5","FC3")
Roi_2 <- c("Fp1","Fp2","AF3","AFz","AF4","F1","Fz","F2","FC1","FC2")
Roi_3 <- c("AF8","F4","F6","F8","FC4","FC6","FT8")
Roi_4 <- c("T7","TP7","C5","C3","CP5","CP3","P7","P5","P3")
Roi_5 <- c("C1","Cz","C2","CP1","CPz","CP2","P1","Pz","P2")
Roi_6 <- c("C4","C6","T8","CP4","CP6","TP8","P4","P6","P8")
Roi_7 <- c("PO9","PO7","PO3","POz","PO4","PO8","PO10","O1","Oz","O2")
Roi_list <- list(Roi_1, Roi_2, Roi_3, Roi_4, Roi_5, Roi_6, Roi_7)

################################
################################ Hyperparameters for the FreSpeD algorithm
f_sampling <- 250 
duration <- 1 # duration of a local time block used for local spectral estimation

windowLen <- duration*f_sampling # length of the window for the local spectral estimation
M_welch <- 5 #floor(0.1*windowLen/2) # Welch regularizer parameter. 

deltaPer <- 0.05 # minimal distance between two change points

#################################
#################################

### Divide eeg recordings of 8 minutes to epochs 

#'@param X: eeg signals (matrix with rows as time and columns as channels)
#'@param epoch_duration: duration of each epoch in seconds
#'@param overlap: overlap between epochs in seconds
divide_to_epochs <- function(X, epoch_duration = 60, overlap = 0, f_sampling = 250) {
  n <- dim(X)[1]  # number of time points in the signal
  
  epoch_samples <- epoch_duration * f_sampling  # Number of samples per epoch
  overlap_samples <- overlap * f_sampling  # Number of samples for overlap
  
  # Calculate the total number of epochs with overlap considered
  n_epochs <- ceiling((n - overlap_samples) / (epoch_samples - overlap_samples))
  
  epochs <- list()  
  for (i in 1:n_epochs) {
    start_idx <- (i - 1) * (epoch_samples - overlap_samples) + 1
    end_idx <- min(start_idx + epoch_samples - 1, n)
    epochs[[i]] <- X[start_idx:end_idx, , drop = FALSE]
  }
  
  return(epochs)
}



get_summary_roi_cp <- function(subject_id = "010002", Data_path = NULL, windowLen, M_welch = 1,f_sampling = 250, deltaPer = 0.05) {
  tryCatch({
    if(!is.null(Data_path)){
      setwd(Data_path)
      open_path <- paste(subject_id ,"_EO.set", sep="")
      closed_path <- paste(subject_id ,"_EC.set", sep="")
    } else {
      open_path <- paste("Data/sub-", subject_id ,"_EO.set", sep="")
      closed_path <- paste("Data/sub-",subject_id ,"_EC.set", sep="")
    }
    # Each file contains informations about the eeg signals and electrodes locations and "???..." 
    eeg_open <- import_set(open_path) 
    eeg_closed <- import_set(closed_path)
    # Get the eeg recordings
    X_o <- eeg_open$signals
    X_c <- eeg_closed$signals
    # Some recordings with closed and open eyes have different electrodes, we keep the common ones
    available_ch <- intersect(colnames(X_o), colnames(X_c))
    # Arnaud summary
     # Determine available frequency bands using your function
    available_bands <- get_available_bands(windowLen = windowLen, f_sampling = f_sampling, M_welch = M_welch)
    
    sum_subject <- tibble(subject = character(),
                       region = integer(),
                       epoch = integer(),
                       eyeState = integer(),
                       cp_sum = list())
    
    i <- 1
    total_steps <- 2*8*7 # Total iterations
    pb <- txtProgressBar(min = 0, max = total_steps, style = 3)  # Create progress bar

    X_o_epochs <- divide_to_epochs(X_o)
    
    for (epoch in 1:length(X_o_epochs)){
      for (region in 1:length(Roi_list)){
        X_tmp <- X_o_epochs[[epoch]][, intersect(Roi_list[[region]],available_ch)]
        cp <- main_function(X_tmp, windowLen = windowLen, M_welch = M_welch, deltaPer = deltaPer, bands_analysis=TRUE, f_sampling = 250)
        
        sum_subject[i,] <- list(subject_id , region, epoch, 1, list(cp))

        setTxtProgressBar(pb, i)  # Update progress bar
        i <- i + 1
      }
    }
    
    X_c_epochs <- divide_to_epochs(X_c)
    
    for (epoch in 1:length(X_c_epochs)){
      for (region in 1:length(Roi_list)){
        X_tmp <- X_c_epochs[[epoch]][, intersect(Roi_list[[region]],available_ch)]
        cp <- main_function(X_tmp, windowLen = windowLen, M_welch = M_welch, deltaPer = deltaPer, bands_analysis=TRUE, f_sampling = 250)

        sum_subject[i,] <- list(subject_id, region, epoch, 0, list(cp))

        setTxtProgressBar(pb, i)  # Update progress bar
        i <- i + 1
      }
    }
    return(sum_subject)
  }, error = function(e) {
    print(paste("Error in get_summary_roi_cp:", e$message)) # Print error message for debugging
    return(NULL)
  })
}


######################################### Arnaud summary 
sub_ids <- read.csv("Data/subject_ids.csv", colClasses = "character")$Subject_ID

##########
# Create an empty list to store results
for (windowLen in c(125,250,500)){
  
  if (windowLen == 125){M_welch = 1}
  if (windowLen == 250){M_welch = 2}
  if (windowLen == 500){M_welch = 4}
  
  sum_all_list <- list()
  eeg_err_id <- c()
  
  for (subject_id in sub_ids) {
    tmp_sum <- get_summary_roi_cp(subject_id,windowLen = windowLen, 
                                  M_welch = M_welch, deltaPer = deltaPer,
                                  f_sampling = f_sampling)
    if (is.null(tmp_sum)) {
      eeg_err_id <- c(eeg_err_id, subject_id)
    } else {
      sum_all_list[[length(sum_all_list) + 1]] <- tmp_sum  # Append to list
    }
    print(paste(subject_id, "done"))
  }
  
  # Combine all results into a tibble
  sum_all <- bind_rows(sum_all_list)
  
  write_rds(sum_all, paste("Output/windowLen",windowLen,"-Welch",M_welch ,".rds",sep=""))
  #sum_subject <- read_rds("???.rds")
}
