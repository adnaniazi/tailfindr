#' Title
#'
#' @param file_path
#' @param save_plots
#' @param show_plots
#' @param save_dir
#' @param plotting_library
#' @param plot_debug
#'
#' @return
#' @export
#'
#' @examples
find_dna_tail_per_read <- function(file_path,
                                   save_plots=FALSE,
                                   show_plots=FALSE,
                                   plot_debug=FALSE,
                                   save_dir='~',
                                   plotting_library='rbokeh'){
    do_plots <- ifelse(save_plots | show_plots, TRUE, FALSE)

    # first read the data and find the tailtype
    data_list <- dna_tailtype_finder(file_path, do_plots)

    read_data <- data_list$read_data
    read_type <- data_list$read_type
    tail_is_valid <- data_list$tail_is_valid
    polya_end <- data_list$polya_end
    polyt_start <- data_list$polyt_start
    samples_per_nt <- read_data$samples_per_nt

    if (!tail_is_valid) {
        return(list(read_id = read_data$read_id,
                    read_type = read_type,
                    tail_is_valid = tail_is_valid,
                    tail_start = NA,
                    tail_end = NA,
                    samples_per_nt = samples_per_nt,
                    tail_length = NA,
                    file_path = file_path))
    }

    # Empirical parameters
    POLY_T_CNDA_THRESHOLD <- 0.20
    POLY_A_CNDA_THRESHOLD <- 0.31
    SPIKE_THRESHOLD <- 2.0
    MOVING_WINDOW_SIZE <- 120
    MAX_GAP_BETWEEN_TAILS <- 1200
    SEC_TAIL_MIN_SIZE <- read_data$samples_per_nt * 15
    SLOPE_THRESHOLD <- 0.20

    # Z-normalize the data
    norm_data <- z_normalize(read_data$raw_data)
    read_length <- length(norm_data)

    # Reverse polyA reads
    if (read_type=='polyA') {
        norm_data <- rev(norm_data)
        polya_start <- read_length - polya_end
        tail_start <- polya_start
    } else {
        tail_start <- polyt_start
    }


    # Recitfy the data (flip anything that is below zero)
    rectified_data <- rectify(norm_data)

    # Windsorize the data (clip anything above a threshold)
    truncated_data <- truncate_spikes(rectified_data,
                                      spike_threshold=SPIKE_THRESHOLD)

    # smoothen the data
    smoothed_data_1 <- left_to_right_sliding_window_cdna_polyt('mean', truncated_data, MOVING_WINDOW_SIZE, 1)
    smoothed_data_2 <- right_to_left_sliding_window_cdna_polyt('mean', truncated_data, MOVING_WINDOW_SIZE, 1)
    smoothed_data_3 <- pmin(smoothed_data_1, smoothed_data_2)
    smoothed_data <- smoothed_data_3

    ### 1. find the slope between consecutive means with a small smoothing window
    window_size <- 10
    step_size <-  10

    start <- tail_start
    end <- min(start+3000, length(norm_data))
    data <- truncated_data[start:end]
    total <- length(data)
    spots <- seq(from = 1, to = (total - window_size), by = step_size)
    mean_data <- vector(mode="numeric", length = length(spots))
    slope <- vector(mode="numeric", length = length(spots)-1)

    for (i in 1:(length(spots)-1)) {
        if (i < length(spots)-2) {
            mean_data[i] <- mean(data[spots[i]:(spots[i+1] - 1)])
            mean_data[i+1] <- mean(data[spots[i+1]:(spots[i+2] - 1)])
            slope[i] = (mean_data[i+1] - mean_data[i]) #/ (spots[i+1] - spots[i])
        } else {
            mean_data[i] <- mean(data[spots[i]:(spots[i+1] - 1)])
            mean_data[i+1] <- mean(data[spots[i+1]:length(data)])
            slope[i] = (mean_data[i+1] - mean_data[i]) #/ (length(data) - spots[i])
        }
    }

    # Refine the tail start
    k <- 1
    precise_tail_start <- NA
    while (k < 20) {
        if ((slope[k] < SLOPE_THRESHOLD) & (slope[k] > -SLOPE_THRESHOLD) &
            (mean_data[k] < SLOPE_THRESHOLD) & (mean_data[k] > -SLOPE_THRESHOLD)) {
            precise_tail_start <- (k-1)*window_size + tail_start
            break
        }
        k <- k + 1
    }

    # Find the tail end
    i <- max(1, k - 1)
    small_glitch_count <- 0
    quit_searching <- FALSE
    tail_end <- NA
    while (i < length(slope)){
        if ((slope[i] < SLOPE_THRESHOLD) & (slope[i] > -SLOPE_THRESHOLD) & (smoothed_data_3[tail_start+i*window_size] < 0.5)) {
            tail_end <- i
            i <- i + 1
        } else {
            j <- i
            while (j < length(slope)) {
                if (j > MAX_GAP_BETWEEN_TAILS/window_size) {
                    quit_searching <- TRUE
                    break
                }
                if ((slope[j] > SLOPE_THRESHOLD) | (slope[j] < -SLOPE_THRESHOLD)) {
                    small_glitch_count <- 0
                    j <- j + 1
                } else {
                    small_glitch_count <- small_glitch_count + 1
                    j <- j + 1
                    if (small_glitch_count > 6){
                        i <- j
                        break
                    }
                }
            }
            if (quit_searching) {
                break
            }
            i <- j
        }
    }

    # if (is.na(tail_end)) {
    #     tail_end <- min(3000 + 50, read_length)
    # } else{
    tail_end <- start + (tail_end+2)*window_size
    # }


    mean_data <- c(rep(NA, times=tail_start),
                   rep(mean_data, each=window_size),
                   rep(NA, times=length(truncated_data)-tail_start-length(mean_data)*window_size))
    slope <- c(rep(NA, times=tail_start),
               rep(slope, each=window_size),
               rep(NA, times=length(truncated_data)-tail_start-length(rep(slope, times=window_size))))


    # reverse the data back to original orientation
    if (do_plots & read_type=='polyA'){
        norm_data <- rev(norm_data)
        truncated_data <- rev(truncated_data)
        smoothed_data_3 <- rev(smoothed_data_3)
        mean_data <- rev(mean_data)
        slope <- rev(slope)
    }

    # correct tail start and ends for polyA reads (that have been reversed)
    if (read_type=='polyA'){
        polya_start <- read_length-tail_end
        polya_end <-  read_length-tail_start
        tail_start <- polya_start
        tail_end <- polya_end
    }

    ##############################################################
    # set the tail_start to precise_tail_start if it was computed
    # N.B.: the location of this code should not be disturbed
    if (!is.na(precise_tail_start)) {
        tail_start <- precise_tail_start
    }
    #############################################################

    # calculate the tail length
    tail_length = (tail_end - tail_start)/read_data$samples_per_nt

    df = data.frame(x=c(1:length(read_data$raw_data)),
                    truncated_data=truncated_data,
                    smoothed_data_3=smoothed_data_3,
                    moves=read_data$moves_sample_wise_vector-3.0,
                    mean_data=mean_data,
                    slope=slope)

    if (!is.na(tail_end)) {
        if (read_type=='polyT') {
            df['poly_t_tail'] <- c(rep(NA, times=tail_start-1),
                                   truncated_data[tail_start:tail_end],
                                   rep(NA, times=(read_length-tail_end)))
        } else {
            df['poly_a_tail'] <- c(rep(NA, times=tail_start-1),
                                   truncated_data[tail_start:tail_end],
                                   rep(NA, times=(read_length-tail_end)))
        }
    }
    filename <- basename(file_path)



    if (plotting_library == 'rbokeh') {
        if (read_type=='polyA'){
            plot_title <- paste('Poly(A) tail  |  ',
                                'Tail length[nt]: ', round(tail_length, 2), '  |  ',
                                'Tail start: ', tail_start, '  |  ',
                                'Tail end: ', tail_end, '  |  ',
                                'Tail duration[Sa]: ', tail_end-tail_start, '  |  ',
                                'Samples per nt: ', samples_per_nt,
                                sep='')
            p <- rbokeh::figure(data=df, width=1000, height=600, legend_location="top_left", title=plot_title)
        } else {
            plot_title <- paste('Poly(T) tail  |  ',
                                'Tail length[nt]: ', round(tail_length, 2), '  |  ',
                                'Tail start: ', tail_start, '  |  ',
                                'Tail end: ', tail_end, '  |  ',
                                'Tail duration[Sa]: ', tail_end-tail_start, '  |  ',
                                'Samples per nt: ', samples_per_nt,
                                sep='')
            p <- rbokeh::figure(data=df, width=1000, height=600, legend_location="top_right", title=plot_title)
        }
        p <- rbokeh::ly_lines(p, truncated_data, width=1.5, color='#4040a1', legend = "Normalized windsorized signal")
        if (plot_debug) {
            p <- rbokeh::ly_lines(p, slope, color='#BF1268', width=2, legend = "Slope of the signal in the search window")
            p <- rbokeh::ly_lines(p, mean_data, color='#F79A14', width=2, legend = "Mean of the signal in the search window")
            p <- rbokeh::ly_lines(p, moves, color='#b2b2b2', width=1, legend = "Moves")
            p <- rbokeh::ly_lines(p, smoothed_data_3, color='#060B54', width=3, legend = "Smoothed signal")
            p <- rbokeh::ly_abline(p, h=SLOPE_THRESHOLD, color = '#50394c', type = 3, width=3, legend = "Slope upper bound")
            p <- rbokeh::ly_abline(p, h=-SLOPE_THRESHOLD, color = '#4040a1',type = 3,  width=3, legend = "Slope lower bound")
            if (!is.na(tail_start)) {
                p <- rbokeh::ly_abline(p, v=tail_start, color = 'orange', type = 2, width=2, legend = "Tail start")
                p <- rbokeh::ly_abline(p, v=tail_end, color = 'red', type = 2, width=2, legend = "Tail end")
            }
        } else {
            if (read_type=='polyT') {
                p <- rbokeh::ly_lines(p, poly_t_tail, color = '#80ced6', legend = "Poly(T) tail")
            } else {
                p <- rbokeh::ly_lines(p, poly_a_tail, color = '#034f84', legend = "Poly(A) tail")
            }
        }
        p <- rbokeh::x_axis(p, label='Sample index')
        p <- rbokeh::y_axis(p, label='z-normalized data values')
        p <- rbokeh::tool_pan(p, dimensions = "width")
        p <- rbokeh::tool_wheel_zoom(p, dimensions = "width")
    }
    if (show_plots) {
        print(p)
    }

    return(list(read_id = read_data$read_id,
                read_type = read_type,
                tail_is_valid = tail_is_valid,
                tail_start = tail_start,
                tail_end = tail_end,
                samples_per_nt = samples_per_nt,
                tail_length = tail_length,
                file_path = file_path))
}
