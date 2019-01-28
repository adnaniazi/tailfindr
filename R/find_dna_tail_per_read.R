#' Title
#'
#' @param file_path
#' @param save_plots
#' @param show_plots
#' @param save_dir
#' @param plotting_library
#' @param plot_debug
#' @param data
#' @param multifast5
#' @param basecalled_with_flipflop
#' @param ...
#' @param read_id_fast5_file
#'
#' @return
#' @export
#'
#' @examples
find_dna_tail_per_read <- function(file_path=NA,
                                   read_id_fast5_file=NA,
                                   data='cdna',
                                   save_plots=F,
                                   show_plots=F,
                                   plot_debug=F,
                                   save_dir=NA,
                                   plotting_library='rbokeh',
                                   multifast5=F,
                                   basecalled_with_flipflop=F,
                                   ...){

    do_plots <- ifelse(save_plots | show_plots, TRUE, FALSE)

    # first read the data and find the tailtype
    if (!multifast5 & !basecalled_with_flipflop) {
        data_list <- dna_tailtype_finder(file_path, plot_debug, data=data)
    } else {
        data_list <- dna_tailtype_finder(plot_debug=plot_debug,
                                         data=data,
                                         multifast5=T,
                                         basecalled_with_flipflop=T,
                                         read_id_fast5_file=read_id_fast5_file)
        file_path <- read_id_fast5_file$fast5_file
    }

    read_data <- data_list$read_data
    event_data <- read_data$event_data
    read_type <- data_list$read_type
    tail_is_valid <- data_list$tail_is_valid
    polya_end <- data_list$polya_end
    polyt_start <- data_list$polyt_start
    samples_per_nt <- read_data$samples_per_nt
    has_precise_boundary <- data_list$has_precise_boundary

    if (!tail_is_valid) {
        return(list(read_id = read_data$read_id,
                    read_type = read_type,
                    tail_is_valid = tail_is_valid,
                    tail_start = NA,
                    tail_end = NA,
                    samples_per_nt = samples_per_nt,
                    tail_length = NA,
                    file_path = file_path,
                    has_precise_boundary = has_precise_boundary))
    }

    # Empirical parameters
    #POLY_T_CNDA_THRESHOLD <- 0.20
    #POLY_A_CNDA_THRESHOLD <- 0.31
    SPIKE_THRESHOLD <- 2.0
    MOVING_WINDOW_SIZE <- 30
    MAX_GAP_BETWEEN_TAILS <- 1200
    #SEC_TAIL_MIN_SIZE <- read_data$samples_per_nt * 15
    SLOPE_THRESHOLD <- 0.20

    # Z-normalize the data
    raw_data <- read_data$raw_data
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

    # Smoothen the data
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

    # Refine the tail start if the boundary we found based on the adaptor is not precise
    k <- 1
    precise_tail_start <- NA
    if (!has_precise_boundary) {
        while (k < 20) {
            if ((slope[k] < SLOPE_THRESHOLD) & (slope[k] > -SLOPE_THRESHOLD) &
                (mean_data[k] < SLOPE_THRESHOLD+0.1) & (mean_data[k] > -SLOPE_THRESHOLD)) {
                precise_tail_start <- (k-1)*window_size + tail_start
                break
            }
            k <- k + 1
        }
    }

    # Find the tail end
    if (!has_precise_boundary) {
        i <- k
    } else {
        i <- 3
        # this is necesarry otherwise the loop below
        # will never find the tail end
    }
    small_glitch_count <- 0
    quit_searching <- FALSE
    tail_end <- NA
    while (i < length(slope)){
        if ((slope[i] < SLOPE_THRESHOLD) & (slope[i] > -SLOPE_THRESHOLD) & (smoothed_data_3[tail_start+i*window_size] < SLOPE_THRESHOLD+0.1)) {
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
                    if (small_glitch_count > 6){ # previously set to 6
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

    # if tail end is not found for whatever reason
    # perhaps due to wrong alignment location of end primer
    # then return
    if (is.na(tail_end)) {
        return(list(read_id = read_data$read_id,
                    read_type = read_type,
                    tail_is_valid = tail_is_valid,
                    tail_start = NA,
                    tail_end = NA,
                    samples_per_nt = samples_per_nt,
                    tail_length = NA,
                    file_path = file_path,
                    has_precise_boundary = has_precise_boundary))
    }

    if (has_precise_boundary){
        tail_end <- start + (tail_end)*window_size
    } else {
        tail_end <- start + (tail_end+4)*window_size
    }

    mean_data <- c(rep(NA, times=tail_start),
                   rep(mean_data, each=window_size),
                   rep(NA, times=length(truncated_data)-tail_start-length(mean_data)*window_size))

    slope <- c(rep(NA, times=tail_start),
               rep(slope, each=window_size),
               rep(NA, times=length(truncated_data)-tail_start-length(rep(slope, times=window_size))))


    # reverse the data back to its original orientation
    # in the case of polyA reads
    if (do_plots & read_type=='polyA'){
        norm_data <- rev(norm_data)
        truncated_data <- rev(truncated_data)
        smoothed_data_3 <- rev(smoothed_data_3)
        mean_data <- rev(mean_data)
        slope <- rev(slope)
    }

    ##############################################################
    # set the tail_start to precise_tail_start if it was computed
    # N.B.: the location of this code should not be disturbed
    if (!is.na(precise_tail_start)) {
        tail_start <- precise_tail_start
    }
    #############################################################

    # correct tail start and ends for polyA reads (that have been reversed)
    if (read_type=='polyA'){
        polya_start <- read_length-tail_end
        polya_end <-  read_length-tail_start
        tail_start <- polya_start
        tail_end <- polya_end
    }

    # calculate the tail length
    tail_length = (tail_end - tail_start)/read_data$samples_per_nt

    df = data.frame(x=c(1:length(read_data$raw_data)),
                    raw_data=raw_data,
                    truncated_data=truncated_data,
                    smoothed_data_3=smoothed_data_3,
                    moves=read_data$moves_sample_wise_vector-3.0,
                    mean_data=mean_data,
                    slope=slope)
    if (multifast5){
        filename <- read_data$read_id
    } else {
        filename <- basename(file_path)
    }

    if (!plot_debug) {
        if (read_type == 'polyA') {
            df['polya_tail'] <- c(rep(NA, times=tail_start-1),
                                   raw_data[tail_start:tail_end],
                                   rep(NA, times=(read_length-tail_end)))
        } else {
            df['polyt_tail'] <- c(rep(NA, times=tail_start-1),
                                   raw_data[tail_start:tail_end],
                                   rep(NA, times=(read_length-tail_end)))
        }
    }


    if (read_type=='polyA'){
        plot_title <- paste('Poly(A) tail  |  ',
                            'Tail length[nt]: ', round(tail_length, 2), '  |  ',
                            'Tail start: ', tail_start, '  |  ',
                            'Tail end: ', tail_end, '  |  ',
                            'Tail duration[Sa]: ', tail_end-tail_start, '  |  ',
                            'Samples per nt: ', samples_per_nt,
                            sep='')
    } else {
        plot_title <- paste('Poly(T) tail  |  ',
                            'Tail length[nt]: ', round(tail_length, 2), '  |  ',
                            'Tail start: ', tail_start, '  |  ',
                            'Tail end: ', tail_end, '  |  ',
                            'Tail duration[Sa]: ', tail_end-tail_start, '  |  ',
                            'Samples per nt: ', samples_per_nt,
                            sep='')
    }
    if (plotting_library == 'rbokeh') {
        if (read_type=='polyA'){
            p <- rbokeh::figure(data=df, width=1000, height=600, legend_location="top_left", title=plot_title)
        } else {
            p <- rbokeh::figure(data=df, width=1000, height=600, legend_location="top_right", title=plot_title)
        }
        if (plot_debug) {
            p <- rbokeh::ly_lines(p, truncated_data, width=1.5, color='#4040a1', legend = "Normalized windsorized signal")
            p <- rbokeh::ly_lines(p, slope, color='#BF1268', width=2, legend = "Slope of the signal in the search window")
            p <- rbokeh::ly_lines(p, mean_data, color='#F79A14', width=2, legend = "Mean of the signal in the search window")
            p <- rbokeh::ly_lines(p, moves, color='#b2b2b2', width=1, legend = "Moves")
            p <- rbokeh::ly_lines(p, smoothed_data_3, color='#060B54', width=3, legend = "Smoothed signal")
            p <- rbokeh::ly_abline(p, h=SLOPE_THRESHOLD, color = '#00B0DF', type = 3, width=2, legend = "Slope upper bound")
            p <- rbokeh::ly_abline(p, h=-SLOPE_THRESHOLD, color = '#FF94CD', type = 3,  width=2, legend = "Slope lower bound")
            if (!is.na(tail_start)) {
                p <- rbokeh::ly_abline(p, v=tail_start, color = 'orange', type = 2, width=2, legend = "Tail start")
                p <- rbokeh::ly_abline(p, v=tail_end, color = 'red', type = 2, width=2, legend = "Tail end")
            }
            p <- rbokeh::y_axis(p, label='z-normalized data values')
        } else {
            p <- rbokeh::ly_lines(p, raw_data, width=1.5, color='#b2b2b2', legend = "Raw data")
            if (read_type=='polyT') {
                p <- rbokeh::ly_lines(p, polyt_tail, color = '#BF1268', legend = "Poly(T) tail")
            } else {
                p <- rbokeh::ly_lines(p, polya_tail, color = '#BF1268', legend = "Poly(A) tail")
            }
            p <- rbokeh::y_axis(p, label='pA')
        }
        p <- rbokeh::x_axis(p, label='Sample index')
        p <- rbokeh::tool_pan(p, dimensions = "width")
        p <- rbokeh::tool_wheel_zoom(p, dimensions = "width")
    } else { # plotting_library == 'ggplot2'

        if (plot_debug) {
            p <- ggplot2::ggplot(data=df, ggplot2::aes(x = x)) +
                ggplot2::geom_line(ggplot2::aes(y = truncated_data), color = '#4040a1')+
                ggplot2::geom_line(ggplot2::aes(y = moves), color = '#b2b2b2') +
                ggplot2::geom_hline(yintercept = SLOPE_THRESHOLD, color = '#00B0DF', linetype = 'dotted') +
                ggplot2::geom_line(ggplot2::aes(y = slope, color = '#BF1268')) +
                ggplot2::geom_hline(yintercept = -SLOPE_THRESHOLD, color = '#FF94CD', linetype = 'dotted') +
                ggplot2::geom_line(ggplot2::aes(y = smoothed_data_3), color='#060B54') +
                ggplot2::geom_line(ggplot2::aes(y = mean_data, color = '#F79A14')) +
                ggplot2::ylab('z-normalized data values')
        } else {
            p <- ggplot2::ggplot(data=df, ggplot2::aes(x = x)) +
                ggplot2::geom_line(ggplot2::aes(y = raw_data), color = '#8E8E8E') +
                ggplot2::ylab('pA')
        }
        if (!is.na(tail_end)) {
            p <- p + ggplot2::geom_line(ggplot2::aes(y = c(rep(NA, times=tail_start-1),
                                                           raw_data[tail_start:tail_end],
                                                           rep(NA, times=(read_length-tail_end)))), color='#BF1268')
        }
        p <- p + ggplot2::ggtitle(plot_title) +
            ggplot2::xlab('Sample index')
    }

    if (show_plots) {
        print(p)
    }

    if (save_plots) {
        filename_png <- paste(filename,'.png', sep='')
        filename_html <- paste(filename,'.html', sep='')
        save_path_png <- file.path(save_dir, 'plots', filename_png, fsep = .Platform$file.sep)
        save_path_html <- file.path(save_dir, 'plots', filename_html, fsep = .Platform$file.sep)
        if (plotting_library == 'ggplot2') {
            ggplot2::ggsave(save_path_png, plot = p, width = 300, height = 70, units = 'mm')
        }
        else {
            #rbokeh::widget2png(p, file = save_path_png)
            rbokeh::rbokeh2html(p, file = save_path_html)
        }
    }

    return(list(read_id = read_data$read_id,
                read_type = read_type,
                tail_is_valid = tail_is_valid,
                tail_start = tail_start,
                tail_end = tail_end,
                samples_per_nt = samples_per_nt,
                tail_length = tail_length,
                file_path = file_path,
                has_precise_boundary = has_precise_boundary))
}


