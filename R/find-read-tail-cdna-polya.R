#' Find Poly(A) tail in a single cDNA read. The read must be a poly(A) read.
#'
#' @param file_path Path of the FAST5 file
#'
#' @return A list of Fast5 file data
#'
#' @examples
#' find_read_tail_cdna_polya('path/to/fast5/file')
find_read_tail_cdna_polya <- function(file_path,
                                      save_plots=FALSE,
                                      show_plots=FALSE,
                                      save_dir='~'){

    # Empirical parameters
    POLY_A_CNDA_THRESHOLD <- 0.31
    POLY_A_CNDA_SPIKE_THRESHOLD <- 2.0
    POLY_A_CDNA_SEC_POLY_A_MAX_GAP <- 1200
    POLY_A_CDNA_MOVING_WINDOW_SIZE <- 120

    # read the FAST5 data
    read_data <- extract_read_data_hdf5r(file_path)
    sampling_rate <- read_data$sampling_rate

    # Z-normalize the data
    norm_data <- z_normalize(read_data$raw_data)

    # Recitfy the data (flip anything that is below zero)
    rectified_data <- rectify(norm_data)

    # Windsorize the data (clip anything above a threshold)
    truncated_data <- truncate_spikes(rectified_data,
                                      spike_threshold=POLY_A_CNDA_SPIKE_THRESHOLD)

    # smoothen the data
    smoothed_data_1 <- left_to_right_sliding_window_cdna_polya('mean', truncated_data, POLY_A_CDNA_MOVING_WINDOW_SIZE, 1)
    smoothed_data_2 <- right_to_left_sliding_window_cdna_polya('mean', truncated_data, POLY_A_CDNA_MOVING_WINDOW_SIZE, 1)
    smoothed_data_3 <- pmin(smoothed_data_1, smoothed_data_2)
    smoothed_data <- smoothed_data_3

    # find intersections with the threshold
    intersections <- smoothed_data < POLY_A_CNDA_THRESHOLD
    rle_intersections <- rle(intersections)

    # merge small intervals in RLE into the bigger intervals
    rle_intersections <- merge_small_rle_intervals(rle_intersections, threshold=20)

    # smoothen the rle, ie remove kinks in it prduced due to the two overlapping windows
    smoothed_rle_intersections <- smoothen_rle_intersections(rle_intersections)
    smoothed_rle_intersections <- smoothen_rle_intersections(smoothed_rle_intersections)

    rle_lengths <- smoothed_rle_intersections$lengths
    rle_values <- smoothed_rle_intersections$values
    rle_indices <- cumsum(rle_lengths)
    len_rle <- length(rle_values)
    read_length <- length(read_data$raw_data)

    # find ploy(A) tail
    cdna_poly_a_read_type <- ''

    # case X
    # file: 1.fast5
    # Run Length Encoding
    # lengths: int [1:2] 946 14784
    # values : logi [1:2] TRUE FALSE
    if (len_rle <= 2){
        cdna_poly_a_read_type <- 'adaptor0'
    }

    # case X
    # file: 5.fast5, 0.fast5
    # from right to left
    # bullshit(FALSE) --> pri_polyA(TRUE) --> bullshit(FALSE) of some length> POLY_A_CDNA_SEC_POLY_A_MAX_GAP
    # $lengths: 1901 9877  734  360
    # $values:TRUE FALSE  TRUE FALSE
    else if (!rle_values[len_rle] && rle_values[(len_rle-1)] &&
             !rle_values[(len_rle-2)]){
        pri_poly_a_start <- rle_indices[(len_rle-2)]
        pri_poly_a_end <- rle_indices[(len_rle-1)]
        cdna_poly_a_read_type <- '..010'

        # find the first secondary tail
        if (len_rle > 4) {
            # if we have a small gap then we a have the first secondary poly-a tail
            # file: 0.fast5
            if ((rle_lengths[(len_rle-2)] < POLY_A_CDNA_SEC_POLY_A_MAX_GAP)&&
                rle_values[(len_rle-3)] && !rle_values[(len_rle-4)]){
                sec1_poly_a_start <- rle_indices[(len_rle-4)]
                sec1_poly_a_end <- rle_indices[(len_rle-3)]
                gap1_start <- sec1_poly_a_end + 1
                gap1_end <- pri_poly_a_start - 1
                cdna_poly_a_read_type <- '..01010'
            }
            # if first secondary tails is found, then find the second secondary tail
            if (len_rle > 6) {
                # if we have a small gap then we a have the first secondary poly-a tail
                # file: 0.fast5
                if ((rle_lengths[(len_rle-4)] < POLY_A_CDNA_SEC_POLY_A_MAX_GAP)&&
                    rle_values[(len_rle-5)] && !rle_values[(len_rle-6)]){
                    sec2_poly_a_start <- rle_indices[(len_rle-6)]
                    sec2_poly_a_end <- rle_indices[(len_rle-5)]
                    gap2_start <- sec2_poly_a_end + 1
                    gap2_end <- sec1_poly_a_start - 1
                    cdna_poly_a_read_type <- '..0101010'
                }
            }
        }
    }

    # case X
    # file: 12.fast5
    # from right to left
    # pri_polyA(TRUE) --> bullshit(FALSE) of some length > POLY_A_CDNA_SEC_POLY_A_MAX_GAP --> adaptor
    # Run Length Encoding
    # lengths: int [1:3] 599 9273 416
    # values : logi [1:3] TRUE FALSE TRUE
    else if (rle_values[len_rle] && !rle_values[(len_rle-1)]
             && rle_values[(len_rle-2)]){
        pri_poly_a_start <- rle_indices[(len_rle-1)]
        pri_poly_a_end <- rle_indices[(len_rle)]
        cdna_poly_a_read_type <- 'adaptor01'
    }

    if (show_plots || save_plots){
        df = data.frame(x=c(1:length(read_data$raw_data)),
                        raw_data=read_data$raw_data,
                        norm_data=norm_data,
                        truncated_data=truncated_data,
                        smoothed_data_1=smoothed_data_1,
                        smoothed_data_2=smoothed_data_2,
                        smoothed_data_3=smoothed_data_3,
                        moves=read_data$moves_sample_wise_vector)

        p <- ggplot2::ggplot(data = df, ggplot2::aes(x = x)) +
             #ggplot2::geom_line(ggplot2::aes(y = norm_data), color='red') +
             ggplot2::geom_line(ggplot2::aes(y = truncated_data), color='blue') +
             ggplot2::geom_line(ggplot2::aes(y = smoothed_data_1), color='black') +
             ggplot2::geom_line(ggplot2::aes(y = smoothed_data_2), color='green') +
             ggplot2::geom_line(ggplot2::aes(y = smoothed_data_3), color='orange') +
             ggplot2::geom_line(ggplot2::aes(y = moves)) +
             ggplot2::geom_hline(yintercept=POLY_A_CNDA_THRESHOLD, color = "black") +
             ggplot2::scale_x_continuous(limits = c(1,
                                                    ceiling(length(smoothed_data_1)/1)))

        if (exists('pri_poly_a_start')) {
            p <- p+ggplot2::geom_line(ggplot2::aes(y = c(rep(NA, times=pri_poly_a_start-1),
                                                         truncated_data[pri_poly_a_start:pri_poly_a_end],
                                                         rep(NA, times=(read_length-pri_poly_a_end)))), color='red')+
            ggplot2::geom_line(ggplot2::aes(y = smoothed_data_3), color='orange')
        }

        if (exists('sec1_poly_a_start')) {
            p <- p+ggplot2::geom_line(ggplot2::aes(y = c(rep(NA, times=sec1_poly_a_start-1),
                                                         truncated_data[sec1_poly_a_start:sec1_poly_a_end],
                                                         rep(NA, times=(read_length-sec1_poly_a_end)))), color='red')+
                ggplot2::geom_line(ggplot2::aes(y = smoothed_data_3), color='orange')
        }

        if (exists('sec2_poly_a_start')) {
            p <- p+ggplot2::geom_line(ggplot2::aes(y = c(rep(NA, times=sec2_poly_a_start-1),
                                                         truncated_data[sec2_poly_a_start:sec2_poly_a_end],
                                                         rep(NA, times=(read_length-sec2_poly_a_end)))), color='red')+
                ggplot2::geom_line(ggplot2::aes(y = smoothed_data_3), color='orange')
        }

        if (show_plots){p}
        if (save_plots){
            filename <- basename(file_path)
            filename <- paste(filename,'.png', sep='')
            dir.create(file.path(save_dir, 'plots', fsep = .Platform$file.sep))
            save_path <- file.path(save_dir, 'plots', filename, fsep = .Platform$file.sep)
            ggplot2::ggsave(save_path, width = 300, height = 70, units = 'mm')
        }

    }

    if (!exists('sec2_poly_a_start')){
        sec2_poly_a_start <- NA
        sec2_poly_a_end <- NA
        gap2_start <- NA
        gap2_end <- NA
    }

    if (!exists('sec1_poly_a_start')){
        sec1_poly_a_start <- NA
        sec1_poly_a_end <- NA
        gap1_start <- NA
        gap1_end <- NA
    }

    if (!exists('pri_poly_a_start')){
        pri_poly_a_start <- NA
        pri_poly_a_end <- NA
    }

    data <- list(pri_poly_a_start=pri_poly_a_start,
                 pri_poly_a_end=pri_poly_a_end,
                 gap1_start=gap1_start,
                 gap1_end=gap1_end,
                 sec1_poly_a_start=sec1_poly_a_start,
                 sec1_poly_a_end=sec1_poly_a_end,
                 gap2_start=gap2_start,
                 gap2_end=gap2_end,
                 sec2_poly_a_start=sec2_poly_a_start,
                 sec2_poly_a_end=sec2_poly_a_end,
                 sampling_rate=sampling_rate,
                 file_path=file_path)

    return(data)
}

