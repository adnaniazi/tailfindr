rm(list = ls())
library(ggplot2)
#source("https://bioconductor.org/biocLite.R")
#biocLite("rhdf5")
library(rhdf5)

# read 54 is problamatic due to a large spike in the leader sequences
# read 20 and 51 has double peaks
# 58 has no polyt but the algortihms still finds a dinstant one TODO

# PRESETS FOR CDNA POLY-T DETECTION
POLY_T_CNDA_THRESHOLD                       <- 0.25  # it is an empirically selected threshold. Anything below this threshold in the smoothed signal is a poly-t candidate
POLY_T_CNDA_ADAPTOR_LENGHT                  <- 300   # unit is in samples (lowest conservative estimate)
POLY_T_CNDA_PRIMARY_SLIDING_WINDOW_SIZE     <- 120   # unit is in samples (24 times the number of raw samples per event). This results is a coarse-windowed signal and is used to detect primary poly-Ts
POLY_T_CNDA_STARTING_GLITCH_SIZE            <- 45    # unit is in samples (There is sharp glitch at the begining of some reads)
POLY_T_CDNA_SECONDARY_POLY_T_MAX_GAP        <- 1000  # if the gap between poly(T)s is more than this number, then they are probably not poly(T)
POLY_T_CNDA_SECONDARY_SLIDING_WINDOW_SIZE   <- 30    # we need a finer window to detect the boundries of shorter secondary poly-Ts. This results in a fine-windowed signal which is used to detect secondary poly-Ts

# -------------------------------------FUNCTIONS--------------------------------

#' z-normalize
#'
#' Apply z-normalization to the time series data
#'
#' @param timeseries Timeseries data that needs to be z-normalized
#'
#' @return z_normalized_timeseries A z-normalized timeseries
#' @examples
#' z_norm(c(1,2,3,4,5))
z_norm <- function(timeseries){
    ts_mean <- mean(timeseries)
    ts_stdev <- sd(timeseries)
    z_normalized_timeseries <- (timeseries - ts_mean)/ts_stdev
    return(z_normalized_timeseries)
}


#' Sliding window
#'
#' Applies a function to window that can slide across a timeseries
#' Acknowledgment: evobiR tool
#'
#' @param FUN Function to apply to each window's data
#' @param data Timeseries data to which to apply the sliding window
#' @param window_size Size of the window
#' @param step_size Step-size of the sliding window
#'
#' @return A signal that has been smoothed using the sliding window. N.B. The
#'   singal is now shorter in lenght
#' @export
#'
#' @examples
#' sliding_window('mean, c(1,2,3,4,5,6,7,8), 2, 1)
#'
sliding_window <- function (FUN, data, window_size, step_size)
{
    total <- length(data)
    spots <- seq(from = 1, to = (total - window_size), by = step_size)
    result <- vector(length = length(spots))
    for (i in 1:length(spots)) {
        result[i] <- match.fun(FUN)(data[spots[i]:(spots[i] + window_size - 1)])
    }
    return(result)
}


#' Remove high-voltage spikes
#'
#' It has been observed the read begins with a high voltage spike that screws up
#' the algorithm. The spike usually within the first 30 samples at the start of
#' a read. This function replaces these spikes with zero.
#'
#' @param data Data which has been z-normalized and then absolutized
#' @param num_samples_to_check The number of starting samples to check for
#'   spikes
#'
#' @return Same data as input sans spikes at the start
#' @export
#'
#' @examples
#' remove_high_voltage_glitch_at_start(data, 30)
#'
remove_high_voltage_glitch_at_start <- function(data, num_samples_to_check){
    # take a portion from the start of the read and remove high-voltage gltichs in it
    data_start <-data[1:num_samples_to_check]
    # find high-voltage spike samples and replace them with 0
    data_start[data_start > 1] = 0
    # put the glitch-free signal back into the raw data and return it
    data[1:num_samples_to_check] <- data_start
    return(data)
}


#' Truncate spikes
#'
#' It has been observed that the signal contains a lot of high-voltage spikes
#' distributed throughout the lenght of the read. This function truncates/clips
#' these spikes to a specified threshold value
#'
#' @param data Data that needs its spikes to be truncated
#' @param spike_threshold Threshold above which spikes will be clipped
#'
#' @return Same data as input but with spikes clipped
#' @export
#'
#' @examples
#' truncate_spikes(data, 2)
#'
truncate_spikes <- function (data, spike_threshold){
    data[data > spike_threshold] = spike_threshold
    return(data)
}


#' Extract FASTA data
#'
#' Extracts read data from basecalled FASTA files. The function can't work on
#' raw reads that haven't been basecalled by Albacore.
#'
#' @param read_path Full path of the FASTA read
#'
#' @return A list of relevant data extracted from the FASTA file
#' @export
#'
#' @examples
extract_read_data <- function(read_path){
    # extract raw data
    f5_tree <- h5ls(read_path)
    raw_reads <- f5_tree[(which(f5_tree == "/Raw/Reads") + 1), 1]
    raw_data <- h5read(read_path, raw_reads)$Signal
    ra <- h5readAttributes(read_path, raw_reads)
    ugk <- h5readAttributes(read_path, 'UniqueGlobalKey/channel_id/')
    
    read_id <- ra$read_id
    read_number <- ra$read_number
    start_time <- ra$start_time
    
    channel_number <- ugk$channel_number
    sampling_rate <- ugk$sampling_rate
    
    # extract events data
    tmpPath <- h5ls(read_path)[which(h5ls(read_path) == "/Analyses/Basecall_1D_000/BaseCalled_template")[1], 1]
    event_data <- h5read(read_path, tmpPath)$Events
    
    # make a vector of moves interpolated for every sample i.e., make a sample-wise or per-sample vector of moves
    if (event_data$start[1] !=0) {
        moves_sample_wise_vector <- c(rep(NA, event_data$start[1]-1),
                                      rep(event_data$move*0.25+1.5, each=event_data$length[1]),
                                      rep( NA, length(raw_data) - (tail(event_data$start, n=1)+event_data$length[1]-1)))
    } else {
        moves_sample_wise_vector <- c(rep(event_data$move*0.25+1.5, each=event_data$length[1]),
                                      rep( NA, length(raw_data) - (tail(event_data$start, n=1)+event_data$length[1])))
        
    }
    
    
    read_data = list(raw_data = raw_data,
                     moves_sample_wise_vector = moves_sample_wise_vector,
                     fast5_event_length = event_data$length[1],
                     read_id = read_id,
                     read_number = read_number,
                     channel_number = channel_number,
                     start_time = start_time,
                     sampling_rate = sampling_rate
    )
    H5close()
    return(read_data)
}

#' Polish Poly-T ends
#'
#' Sometimes the algorithms fails to capture the full-length of poly-T towards
#' the right end due to shift in the averaged signal. In that case, we extend
#' the Poly-T tail using the no-change in moves criteria. Currently we do this
#' only for primary poly-Ts
#'
#' @param primary_poly_t_coords A vector of poly-T coordinates for prmimary
#'   poly-T
#' @param moves_sample_wise_vector A vector of per sample moves
#'
#' @return A vector of poly-T coordinates with the right end of the poly-T
#'   polished with moves
#' @export
#'
#' @examples
#'
polish_end_of_poly_T <- function(primary_poly_t_coords, moves_sample_wise_vector){
    i <- primary_poly_t_coords[2]
    while (TRUE) {
        if (moves_sample_wise_vector[i] == 1.5){
            i <- i + 1
        }
        else {
            break
        }
    }
    primary_poly_t_coords[2] <- i
    return(primary_poly_t_coords)
}


#' Polish the start of poly-T tail
#'
#' Sometimes the algorithms fails to capture the full-length of poly-T towards
#' the left end due to shift in the averaged signal. In that case, we extend the
#' Poly-T tail toward the left side using the no-change in moves criteria.
#' Currently we do this only for primary poly-Ts
#'
#' @param primary_poly_t_coords A vector of poly-T coordinates for prmimary
#'   poly-T
#' @param moves_sample_wise_vector A vector of per sample moves
#'
#' @return A vector of poly-T coordinates with the left end of the poly-T
#'   polished with moves
#' @export
#'
#' @examples
polish_start_of_poly_T <- function(primary_poly_t_coords, moves_sample_wise_vector){
    i <- primary_poly_t_coords[1]
    while (TRUE) {
        if (moves_sample_wise_vector[i] == 1.5){
            i <- i - 1
        }
        else {
            break
        }
    }
    primary_poly_t_coords[1] <- i
    return(primary_poly_t_coords)
}


#' Find total moves in a region
#'
#' This function finds the total number of moves in a range.
#'
#' @param moves_sample_wise_vector A per sample vector of moves
#' @param coordinates A vector of range in which to find the total number of moves
#' @param samples_per_fast5_event Number of samples per event in FAST5 file
#'
#' @return
#' @export
#'
#' @examples
find_moves_in_a_given_length_of_samples <- function(moves_sample_wise_vector, coordinates, samples_per_fast5_event){
    moves_sample_wise_vector <- (moves_sample_wise_vector - 1.5) * 4 # transforming from 1.5, 1.75, and 2 back to 0, 1, and 2
    interval_length <- length(coordinates[1]:coordinates[2])
    moves <- sum(moves_sample_wise_vector[coordinates[1]:coordinates[2]], na.rm = TRUE)/samples_per_fast5_event
    data = list(coordinates <- coordinates,
                interval_length <- interval_length,
                moves_in_the_interval <- moves)
    return(data)
}


remove_spikes_in_rle_intersections <- function(intersections_rle_coarse, merge_anything_below_this_threshold){
    l <- intersections_rle_coarse$lengths
    v <- intersections_rle_coarse$values
    
    less_than <- intersections_rle_coarse$lengths < merge_anything_below_this_threshold
    
    less_than_threshold_idx <- which(intersections_rle_coarse$lengths %in% intersections_rle_coarse$lengths[less_than])
    
    for (i in less_than_threshold_idx){
        if (i>1){
            l[(i-1)] <- l[(i-1)]+l[i]
        }
    }
    l_backup <- l
    l <- l[which(! l %in% l[less_than_threshold_idx])]
    v <- v[which(! l_backup %in% l_backup[less_than_threshold_idx])]
    intersections_rle_coarse = list(lengths=l,
                                    values=v)
    return(intersections_rle_coarse)
}

main_find_poly_t <- function(path2file, image_save_path='~', save_image=FALSE, display_image=FALSE){
    # ---------------------------------------CODE-----------------------------------
    #path2file <- '/Users/adnaniazi/Documents/phd/code/ploy-A-browser/cdna_polyt_fast5s/0.fast5'
    #path2file <- '/Users/adnaniazi/Dropbox/r_code/cdna_polyt_fast5s/1.fast5'
    read_data <- extract_read_data(path2file)
    raw_data <- read_data$raw_data
    moves_sample_wise_vector <- read_data$moves_sample_wise_vector
    raw_data_znorm <- z_norm(raw_data)
    raw_data_znorm_abs <-  abs(raw_data_znorm)
    processed_raw_data <- remove_high_voltage_glitch_at_start(raw_data_znorm_abs,
                                                              POLY_T_CNDA_STARTING_GLITCH_SIZE)
    
    # anything above 2 standard deviation is abnormal. Truncate these samples
    processed_raw_data <- truncate_spikes(processed_raw_data, 2)
    
    # make a coarse and fine-windowed signal to detect the primary and secondary poly-Ts
    # coarse-grained sliding window to find the primary poly-T
    coarse_sliding_windowed_signal <- sliding_window('mean',
                                                     processed_raw_data,
                                                     120, 1)
    # fine-grained sliding window to find the secondary poly-T
    fine_sliding_windowed_signal <- sliding_window('mean',
                                                   processed_raw_data,
                                                   30, 1)
    
    # Find intersection of the smoothed windowed curves with the POLY-T threshold.
    # Anything below the threshold is a candiate poly-T
    intersections_coarse <-  coarse_sliding_windowed_signal < POLY_T_CNDA_THRESHOLD
    intersections_fine <-  fine_sliding_windowed_signal < POLY_T_CNDA_THRESHOLD
    
    intersections_rle_coarse <- rle(intersections_coarse)
    intersections_rle_fine <- rle(intersections_fine)
    intersections_rle_coarse <- remove_spikes_in_rle_intersections(intersections_rle_coarse, 10)
    
    # Terminologies:
    # low-level == anything below the threshold (these are poly-T candiates)
    # high-level == anything above the threshold (these are not poly-Ts)
    
    # CASE 1:
    # Example Fast5 file: 8.fast5
    # - starts with low-level poly-t looking leader
    # - followed by high-level adaptor
    # - followed by the actual low-level poly-t istself
    # - followed by high-level cdna sequence
    primary_poly_t_coords <- c(NA, NA)
    if (identical(intersections_rle_coarse$values[1:4], c(TRUE, FALSE, TRUE, FALSE)) &&
        intersections_rle_coarse$lengths[2] > POLY_T_CNDA_ADAPTOR_LENGHT) {
        
        poly_t_type <- 1
        
        primary_poly_t_coords <- c(sum(intersections_rle_coarse$lengths[1:2], -POLY_T_CNDA_PRIMARY_SLIDING_WINDOW_SIZE/4),
                                   sum(intersections_rle_coarse$lengths[1:3], POLY_T_CNDA_PRIMARY_SLIDING_WINDOW_SIZE/2)) + #  POLY_T_CNDA_PRIMARY_SLIDING_WINDOW_SIZE/2) expands the poly-T
            POLY_T_CNDA_PRIMARY_SLIDING_WINDOW_SIZE/4 # it shifts the poly-T towords the right to correct location
        
        # CORNER CASE: Example file: 50.fast5
        # Sometime the poly-t detected by the above algorithm starts at the correct location
        # but ends prematurely. Use the moves to refine the end of the primary poly-ts
        # 0 move = 1.5, 1 move = 1.75, and 2 moves = 2 in moves_sample_wise_vector
        primary_poly_t_coords <- polish_end_of_poly_T(primary_poly_t_coords, moves_sample_wise_vector)
        primary_poly_t_coords <- polish_start_of_poly_T(primary_poly_t_coords, moves_sample_wise_vector)
        poly_t_adaptor_coords <- c(intersections_rle_coarse$lengths[1][1], primary_poly_t_coords[1]-1)
        poly_t_adaptor_data <- find_moves_in_a_given_length_of_samples(moves_sample_wise_vector, poly_t_adaptor_coords, read_data$fast5_event_length)
        primary_poly_t_data <- find_moves_in_a_given_length_of_samples(moves_sample_wise_vector, primary_poly_t_coords, read_data$fast5_event_length)
        
        # FIND FIRST SECONDARY POLY-T
        # Example file: 5.fast5
        #check if there is gap and that this gap is less than the maximum allowed gap between primary and secondary poly-T
        if (!intersections_rle_coarse$values[4] && intersections_rle_coarse$values[5] && intersections_rle_coarse$lengths[4]<POLY_T_CDNA_SECONDARY_POLY_T_MAX_GAP) {
            
            #coarse-grained boundries of first gap and first secondary poly-T
            first_gap_coords_coarse <- c(primary_poly_t_coords[2]+1, primary_poly_t_coords[2] + 1 + intersections_rle_coarse$lengths[4])
            first_secondary_poly_t_coords_coarse <- c(first_gap_coords_coarse[2]+1, first_gap_coords_coarse[2] + 1 + intersections_rle_coarse$lengths[5])
            
            #fine-grained boundries of first gap and first secondary poly-T
            intersections_rle_fine_cumsum <- cumsum(intersections_rle_fine$lengths)
            
            gap1_end <- tail(intersections_rle_fine_cumsum[intersections_rle_fine_cumsum < first_gap_coords_coarse[2]], n=1)
            first_gap_coords_fine <- c(primary_poly_t_coords[2]+1, gap1_end)
            
            sec_polyt1_end <- tail(intersections_rle_fine_cumsum[intersections_rle_fine_cumsum < first_secondary_poly_t_coords_coarse[2]], n=1)
            first_secondary_poly_t_coords_fine <- c(first_gap_coords_fine[2]+1, sec_polyt1_end+POLY_T_CNDA_SECONDARY_SLIDING_WINDOW_SIZE/6)
            first_secondary_poly_t_coords_fine <- polish_end_of_poly_T(first_secondary_poly_t_coords_fine, moves_sample_wise_vector)
            
            first_secondary_poly_t_data <- find_moves_in_a_given_length_of_samples(moves_sample_wise_vector, first_secondary_poly_t_coords_fine, read_data$fast5_event_length)
            first_gap_coords_data <- find_moves_in_a_given_length_of_samples(moves_sample_wise_vector, first_gap_coords_fine, read_data$fast5_event_length)
            
            # If YOU FIND FIRST SECONDARY POLY_T, THEN LOOK FOR A SCEOND ONE TOO
            # Example file: 57.fast5
            # check if there is gap and that this gap is less than the maximum allowed gap between primary and secondary poly-T
            if (intersections_rle_coarse$values[5] && !intersections_rle_coarse$values[6] && intersections_rle_coarse$lengths[6]<POLY_T_CDNA_SECONDARY_POLY_T_MAX_GAP){
                
                # coarse-grained boundries of second gap and second secondary poly-T
                second_gap_coords_coarse <- c(first_secondary_poly_t_coords_coarse[2]+1, first_secondary_poly_t_coords_coarse[2] + 1 + intersections_rle_coarse$lengths[5])
                second_secondary_poly_t_coords_coarse <- c(second_gap_coords_coarse[2]+1, second_gap_coords_coarse[2] + 1 + intersections_rle_coarse$lengths[6])
                
                #fine-grained boundries of first gap and first secondary poly-T
                intersections_rle_fine_cumsum <- cumsum(intersections_rle_fine$lengths)
                
                # refine the end of the second gap by setting to the intersection of fine-windowed signal
                gap2_end_left <- tail(intersections_rle_fine_cumsum[intersections_rle_fine_cumsum < second_gap_coords_coarse[2]], n=1)
                gap2_end_right <- head(intersections_rle_fine_cumsum[intersections_rle_fine_cumsum > second_gap_coords_coarse[2]], n=1)
                if ( abs(gap2_end_left - second_gap_coords_coarse[2]) < abs(gap2_end_right - second_gap_coords_coarse[2]) ){
                    second_gap_coords_fine <- c(first_secondary_poly_t_coords_fine[2]+1, gap2_end_left)
                } else {
                    second_gap_coords_fine <- c(first_secondary_poly_t_coords_fine[2]+1, gap2_end_right)
                }
                
                sec_polyt2_end <- tail(intersections_rle_fine_cumsum[intersections_rle_fine_cumsum < second_secondary_poly_t_coords_coarse[2]], n=1)
                second_secondary_poly_t_coords_fine <- c(second_gap_coords_fine[2]+1, sec_polyt2_end+POLY_T_CNDA_SECONDARY_SLIDING_WINDOW_SIZE/6)
                second_secondary_poly_t_coords_fine <- polish_end_of_poly_T(second_secondary_poly_t_coords_fine, moves_sample_wise_vector)
                
                second_secondary_poly_t_data <- find_moves_in_a_given_length_of_samples(moves_sample_wise_vector, second_secondary_poly_t_coords_fine, read_data$fast5_event_length)
                second_gap_coords_data <- find_moves_in_a_given_length_of_samples(moves_sample_wise_vector, second_gap_coords_fine, read_data$fast5_event_length)
            }
        }
        
        # CASE 2:
        # Example Fast5 file: 1.fast5
        # - starts with a high-level adaptor
        # - followed by the actual low-level poly-t istself
        # - followed by high-level cdna sequence
    }else if (identical(intersections_rle_coarse$values[1:3], c(FALSE, TRUE, FALSE)) &&
              intersections_rle_coarse$lengths[1] > POLY_T_CNDA_ADAPTOR_LENGHT) {
        
        poly_t_type <- 2
        
        primary_poly_t_coords <- c(sum(intersections_rle_coarse$lengths[1:1], -POLY_T_CNDA_PRIMARY_SLIDING_WINDOW_SIZE/4),
                                   sum(intersections_rle_coarse$lengths[1:2], POLY_T_CNDA_PRIMARY_SLIDING_WINDOW_SIZE/2)) + #  POLY_T_CNDA_PRIMARY_SLIDING_WINDOW_SIZE/2) expands the poly-T
            POLY_T_CNDA_PRIMARY_SLIDING_WINDOW_SIZE/4 # it shifts the poly-T towords the right to correct location
        
        # CORNER CASE: Example file: 50.fast5
        # Sometime the poly-t detected by the above algorithm starts at the correct location
        # but ends prematurely. Use the moves to refine the end of the primary poly-ts
        # 0 move = 1.5, 1 move = 1.75, and 2 moves = 2 in moves_sample_wise_vector
        primary_poly_t_coords <- polish_end_of_poly_T(primary_poly_t_coords, moves_sample_wise_vector)
        primary_poly_t_coords <- polish_start_of_poly_T(primary_poly_t_coords, moves_sample_wise_vector)
        poly_t_adaptor_coords <- c(0, primary_poly_t_coords[1]-1)
        poly_t_adaptor_data <- find_moves_in_a_given_length_of_samples(moves_sample_wise_vector, poly_t_adaptor_coords, read_data$fast5_event_length)
        primary_poly_t_data <- find_moves_in_a_given_length_of_samples(moves_sample_wise_vector, primary_poly_t_coords, read_data$fast5_event_length)
        
        # FIND THE FIRST SECONDARY POLY-Ts
        # Example file:
        #check if there is gap and that this gap is less than the maximum allowed gap between primary and secondary poly-T
        if (!intersections_rle_coarse$values[3] && intersections_rle_coarse$values[4] && intersections_rle_coarse$lengths[3]<POLY_T_CDNA_SECONDARY_POLY_T_MAX_GAP){
            
            #coarse-grained boundries of first gap and first secondary poly-T
            first_gap_coords_coarse <- c(primary_poly_t_coords[2]+1, primary_poly_t_coords[2] + 1 + intersections_rle_coarse$lengths[3])
            first_secondary_poly_t_coords_coarse <- c(first_gap_coords_coarse[2]+1, first_gap_coords_coarse[2] + 1 + intersections_rle_coarse$lengths[4])
            
            #fine-grained boundries of first gap and first secondary poly-T
            intersections_rle_fine_cumsum <- cumsum(intersections_rle_fine$lengths)
            
            gap1_end <- tail(intersections_rle_fine_cumsum[intersections_rle_fine_cumsum < first_gap_coords_coarse[2]], n=1)
            first_gap_coords_fine <- c(primary_poly_t_coords[2]+1, gap1_end)
            
            sec_polyt1_end <- tail(intersections_rle_fine_cumsum[intersections_rle_fine_cumsum < first_secondary_poly_t_coords_coarse[2]], n=1)
            first_secondary_poly_t_coords_fine <- c(first_gap_coords_fine[2]+1, sec_polyt1_end+POLY_T_CNDA_SECONDARY_SLIDING_WINDOW_SIZE/6)
            first_secondary_poly_t_coords_fine <- polish_end_of_poly_T(first_secondary_poly_t_coords_fine, moves_sample_wise_vector)
            
            first_secondary_poly_t_data <- find_moves_in_a_given_length_of_samples(moves_sample_wise_vector, first_secondary_poly_t_coords_fine, read_data$fast5_event_length)
            first_gap_coords_data <- find_moves_in_a_given_length_of_samples(moves_sample_wise_vector, first_gap_coords_fine, read_data$fast5_event_length)
            
            # If YOU FIND FIRST SECONDARY POLY_T, THEN LOOK FOR A SCEOND ONE TOO
            # Example file:
            # check if there is gap and that this gap is less than the maximum allowed gap between primary and secondary poly-T
            if (intersections_rle_coarse$values[4] && !intersections_rle_coarse$values[5] && intersections_rle_coarse$lengths[5]<POLY_T_CDNA_SECONDARY_POLY_T_MAX_GAP){
                
                # coarse-grained boundries of second gap and second secondary poly-T
                second_gap_coords_coarse <- c(first_secondary_poly_t_coords_coarse[2]+1, first_secondary_poly_t_coords_coarse[2] + 1 + intersections_rle_coarse$lengths[4])
                second_secondary_poly_t_coords_coarse <- c(second_gap_coords_coarse[2]+1, second_gap_coords_coarse[2] + 1 + intersections_rle_coarse$lengths[5])
                
                #fine-grained boundries of first gap and first secondary poly-T
                intersections_rle_fine_cumsum <- cumsum(intersections_rle_fine$lengths)
                
                # refine the end of the second gap by setting to the intersection of fine-windowed signal
                gap2_end_left <- tail(intersections_rle_fine_cumsum[intersections_rle_fine_cumsum < second_gap_coords_coarse[2]], n=1)
                gap2_end_right <- head(intersections_rle_fine_cumsum[intersections_rle_fine_cumsum > second_gap_coords_coarse[2]], n=1)
                if ( abs(gap2_end_left - second_gap_coords_coarse[2]) < abs(gap2_end_right - second_gap_coords_coarse[2]) ){
                    second_gap_coords_fine <- c(first_secondary_poly_t_coords_fine[2]+1, gap2_end_left)
                } else {
                    second_gap_coords_fine <- c(first_secondary_poly_t_coords_fine[2]+1, gap2_end_right)
                }
                
                sec_polyt2_end <- tail(intersections_rle_fine_cumsum[intersections_rle_fine_cumsum < second_secondary_poly_t_coords_coarse[2]], n=1)
                second_secondary_poly_t_coords_fine <- c(second_gap_coords_fine[2]+1, sec_polyt2_end+POLY_T_CNDA_SECONDARY_SLIDING_WINDOW_SIZE/6)
                second_secondary_poly_t_coords_fine <- polish_end_of_poly_T(second_secondary_poly_t_coords_fine, moves_sample_wise_vector)
                
                second_secondary_poly_t_data <- find_moves_in_a_given_length_of_samples(moves_sample_wise_vector, second_secondary_poly_t_coords_fine, read_data$fast5_event_length)
                second_gap_coords_data <- find_moves_in_a_given_length_of_samples(moves_sample_wise_vector, second_gap_coords_fine, read_data$fast5_event_length)
            }
        }
    } else {
        # The read does not conform to the expectation. No poly-T tail can be found in this case
        poly_t_type <- 3
    }
    
    
    if (poly_t_type != 3) {
        # find length and moves of post poly-T region
        if (exists('second_secondary_poly_t_data')){
            post_poly_t_coords <- c(second_secondary_poly_t_data[[1]][2]+1, length(processed_raw_data))
            post_poly_t_data <- find_moves_in_a_given_length_of_samples(moves_sample_wise_vector, post_poly_t_coords, read_data$fast5_event_length)
        } else if (exists('first_secondary_poly_t_data')){
            post_poly_t_coords <- c(first_secondary_poly_t_data[[1]][2]+1, length(processed_raw_data))
            post_poly_t_data <- find_moves_in_a_given_length_of_samples(moves_sample_wise_vector, post_poly_t_coords, read_data$fast5_event_length)
        }  else if (exists('primary_poly_t_data')){
            post_poly_t_coords <- c(primary_poly_t_data[[1]][2]+1, length(processed_raw_data))
            post_poly_t_data <- find_moves_in_a_given_length_of_samples(moves_sample_wise_vector, post_poly_t_coords, read_data$fast5_event_length)
        }
    } else {
        post_poly_t_data <- list(c(NA,NA), NA, NA)
    }
    
    
    if ((poly_t_type != 3) && save_image){
        df = data.frame(x=c(1:length(processed_raw_data)),
                        y1=processed_raw_data,
                        y2=c(coarse_sliding_windowed_signal, rep(NA, length(processed_raw_data)-length(coarse_sliding_windowed_signal))),
                        y3=c(rep(NA, primary_poly_t_coords[1]-1), processed_raw_data[primary_poly_t_coords[1]:primary_poly_t_coords[2]], rep(NA, length(processed_raw_data)-primary_poly_t_coords[2])),
                        moves_sample_wise_vector = moves_sample_wise_vector
        )
        p <- ggplot(data = df, aes(x=x, y=x)) +
            geom_line(aes(y=y1), color='red') +
            geom_line(aes(y=y3), color = 'blue') +
            geom_line(aes(y=moves_sample_wise_vector), color = 'orange') +
            geom_line(aes(y=y2)) +
            geom_hline(yintercept = POLY_T_CNDA_THRESHOLD)
        
        if (exists('first_gap_coords_fine')){
            sec_df = data.frame(x=c(1:length(processed_raw_data)),
                                windowed_data=c(coarse_sliding_windowed_signal, rep(NA, length(processed_raw_data)-length(coarse_sliding_windowed_signal))),
                                first_gap=c(rep(NA, primary_poly_t_coords[2]), processed_raw_data[first_gap_coords_fine[1]:first_gap_coords_fine[2]], rep(NA, length(processed_raw_data)-first_gap_coords_fine[2])),
                                first_secondary_poly_t=c(rep(NA, first_gap_coords_fine[2]), processed_raw_data[first_secondary_poly_t_coords_fine[1]:first_secondary_poly_t_coords_fine[2]], rep(NA, length(processed_raw_data)-first_secondary_poly_t_coords_fine[2])))
            
            
            p <- p + geom_line(data = sec_df, aes(y=first_gap), color = 'green') +
                geom_line(data = sec_df, aes(y=first_secondary_poly_t), color = 'blue') +
                geom_line(data = sec_df, aes(y=windowed_data), color = 'black') +
                geom_line(data = sec_df, aes(y=moves_sample_wise_vector), color = 'orange')
        }
        
        if (exists('second_gap_coords_fine')){
            sec_df = data.frame(x=c(1:length(processed_raw_data)),
                                windowed_data=c(coarse_sliding_windowed_signal, rep(NA, length(processed_raw_data)-length(coarse_sliding_windowed_signal))),
                                second_gap=c(rep(NA, first_secondary_poly_t_coords_fine[2]), processed_raw_data[second_gap_coords_fine[1]:second_gap_coords_fine[2]], rep(NA, length(processed_raw_data)-second_gap_coords_fine[2])),
                                secondary_secondary_poly_t=c(rep(NA, second_gap_coords_fine[2]), processed_raw_data[second_secondary_poly_t_coords_fine[1]:second_secondary_poly_t_coords_fine[2]], rep(NA, length(processed_raw_data)-second_secondary_poly_t_coords_fine[2])))
            
            
            p <- p + geom_line(data = sec_df, aes(y=second_gap), color = 'green') +
                geom_line(data = sec_df, aes(y=secondary_secondary_poly_t), color = 'blue') +
                geom_line(data = sec_df, aes(y=windowed_data), color = 'black')
        }
        
        filename <- basename(path2file)
        filename <- paste(filename,'.png', sep='')
        filepath <- file.path(image_save_path, filename, fsep = .Platform$file.sep)
        ggsave(filepath, width = 300, height = 70, units = 'mm')
    } else {
        p <- -1
    }
    
    if (display_image){
        print(p)
    }
    
    ## RETURNS
    if (!exists('second_secondary_poly_t_data')) {
        second_secondary_poly_t_data <- list(c(NA,NA), NA, NA)
        second_gap_coords_data <- list(c(NA,NA), NA, NA)
    }
    if (!exists('first_secondary_poly_t_data')) {
        first_secondary_poly_t_data <- list(c(NA,NA), NA, NA)
        first_gap_coords_data <- list(c(NA,NA), NA, NA)
    }
    if (!exists('primary_poly_t_data')) {
        primary_poly_t_data <- list(c(NA,NA), NA, NA)
    }
    if (!exists('poly_t_adaptor_data')) {
        poly_t_adaptor_data <- list(c(NA,NA), NA, NA)
    }
    
    
    poly_t_data_list <- list(poly_t_adaptor_data = poly_t_adaptor_data,
                             primary_poly_t_data = primary_poly_t_data,
                             first_gap_coords_data = first_gap_coords_data,
                             first_secondary_poly_t_data = first_secondary_poly_t_data,
                             second_gap_coords_data = second_gap_coords_data,
                             second_secondary_poly_t_data = second_secondary_poly_t_data,
                             post_poly_t_data = post_poly_t_data,
                             fast5_event_length = read_data$fast5_event_length,
                             read_id = read_data$read_id,
                             read_number = read_data$read_number,
                             channel_number = read_data$channel_number,
                             start_time = read_data$start_time,
                             sampling_rate = read_data$sampling_rate,
                             poly_t_type = poly_t_type,
                             read_path = path2file
    )
    
    return(poly_t_data_list)
    
} #END of Main Function

update_dataframe <- function(df, pt){
    
    row <- data.frame(adaptor_start=unlist(pt[[1]][1])[1], adaptor_end=unlist(pt[[1]][1])[2], adaptor_length=unlist(pt[[1]][2]), adaptor_moves=unlist(pt[[1]][3]),
                      pri_poly_t_start=unlist(pt[[2]][1])[1], pri_poly_t_end=unlist(pt[[2]][1])[2], pri_poly_t_length=unlist(pt[[2]][2]), pri_poly_t_moves=unlist(pt[[2]][3]),
                      first_gap_start=unlist(pt[[3]][1])[1], first_gap_end=unlist(pt[[3]][1])[2], first_gap_length=unlist(pt[[3]][2]), first_gap_t_moves=unlist(pt[[3]][3]),
                      first_sec_poly_t_start=unlist(pt[[4]][1])[1], first_sec_poly_t_end=unlist(pt[[4]][1])[2], first_sec_poly_t_length=unlist(pt[[4]][2]), first_sec_poly_t_moves=unlist(pt[[4]][3]),
                      second_gap_start=unlist(pt[[5]][1])[1], second_gap_end=unlist(pt[[5]][1])[2], second_gap_length=unlist(pt[[5]][2]), second_gap_t_moves=unlist(pt[[5]][3]),
                      second_sec_poly_t_start=unlist(pt[[6]][1])[1], second_sec_poly_t_end=unlist(pt[[6]][1])[2], second_sec_poly_t_length=unlist(pt[[6]][2]), second_sec_poly_t_moves=unlist(pt[[6]][3]),
                      post_poly_t_start=unlist(pt[[7]][1])[1], post_poly_t_end=unlist(pt[[7]][1])[2], post_poly_t_length=unlist(pt[[7]][2]), post_poly_t_moves=unlist(pt[[7]][3]),
                      fast5_event_length=unlist(pt[[8]]), read_id=unlist(pt[[9]]), read_number=unlist(pt[[10]]),
                      channel_number=unlist(pt[[11]]), start_time=unlist(pt[[12]]), sampling_rate=unlist(pt[[13]]),
                      poly_t_type=unlist(pt[[14]]), read_path=unlist(pt[[15]]))
    
    poly_t_df <- rbind(df, row)
    return(poly_t_df)
}


main_find_poly_t_wrapper <- function(filepath, image_save_path=save_path, save_image=FALSE, display_image=FALSE){
    output <- tryCatch({
        main_find_poly_t(filepath, image_save_path=save_path, save_image=FALSE, display_image=FALSE)
    },
    error = function(cond){
        return(-9)
    }
    )
}

make_empty_polt_t_data_dataframe <- function(){
    col_names_poly_t_df <- c('adaptor_start', 'adaptor_end', 'adaptor_length', 'adaptor_moves',
                             'pri_poly_t_start', 'pri_poly_t_end', 'pri_poly_t_length', 'pri_poly_t_moves',
                             'first_gap_start', 'first_gap_end', 'first_gap_length', 'first_gap_t_moves',
                             'first_sec_poly_t_start', 'first_sec_poly_t_end', 'first_sec_poly_t_length', 'first_sec_poly_t_moves',
                             'second_gap_start', 'second_gap_end', 'second_gap_length', 'second_gap_t_moves',
                             'second_sec_poly_t_start', 'second_sec_poly_t_end', 'second_sec_poly_t_length', 'second_sec_poly_t_moves',
                             'post_poly_t_start', 'post_poly_t_end', 'post_poly_t_length', 'post_poly_t_moves',
                             'fast5_event_length', 'read_id', 'read_number',
                             'channel_number', 'start_time', 'sampling_rate',
                             'poly_t_type', 'read_path')
    poly_t_df <- data.frame(names(col_names_poly_t_df))
    return(poly_t_df)
}

######################################## CALL MAIN ########################################


fast5_dir <- '/export/valenfs/data/processed_data/MinION/20180516_1429_polya_cdna_shield_run1/demultiplexed_gfp_spikeins/poly_t_reads_100_percent_match_threshold/fished_out_reads/_10bp_reads'
save_path <- '/export/valenfs/data/processed_data/MinION/20180516_1429_polya_cdna_shield_run1/demultiplexed_gfp_spikeins/poly_t_reads_100_percent_match_threshold/fished_out_reads/_10bp_reads_analysis'
fast5_files_list <- list.files(path = fdir, pattern = "\\.fast5$", recursive = TRUE, full.names = TRUE)
poly_t_df <- make_empty_polt_t_data_dataframe()

i <- 1
for (f in fast5_files_list){
    print(f)
    pt <- main_find_poly_t_wrapper(f, image_save_path=save_path, save_image=TRUE, display_image=FALSE)
    if (!is.list(pt)) {
        pt <- list(list(c(NA,NA),NA,NA),list(c(NA,NA),NA,NA),list(c(NA,NA),NA,NA),
                   list(c(NA,NA),NA,NA),list(c(NA,NA),NA,NA),list(c(NA,NA),NA,NA),
                   list(c(NA,NA),NA,NA),NA, NA, NA, NA,NA, NA, -4, f)
    }
    poly_t_df <- update_dataframe(poly_t_df, pt)
    print(i)
    i <- i + 1
}

data.table::fwrite(poly_t_df, "/export/valenfs/data/processed_data/MinION/20180516_1429_polya_cdna_shield_run1/demultiplexed_gfp_spikeins/poly_t_reads_100_percent_match_threshold/fished_out_reads/_10bp_reads_analysis/bp10.csv")

###########################################################################################
