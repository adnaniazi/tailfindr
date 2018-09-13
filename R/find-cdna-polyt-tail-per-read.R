#' Find Poly(T) tail in a single cDNA read. The read must be a poly(A) read.
#'
#' @param file_path Path of the FAST5 file
#'
#' @return A list of Fast5 file data
#' @export
#'
#' @examples
#' find_cdna_polyt_tail_per_read('path/to/fast5/file')
find_cdna_polyt_tail_per_read <- function(file_path,
                                          save_plots=FALSE,
                                          show_plots=FALSE,
                                          save_dir='~'){

    # Empirical parameters
    POLY_T_CNDA_THRESHOLD <- 0.31
    POLY_T_CNDA_SPIKE_THRESHOLD <- 2.0
    POLY_T_CDNA_SEC_POLY_A_MAX_GAP <- 1200
    POLY_T_CDNA_MOVING_WINDOW_SIZE <- 120

    # read the FAST5 data
    read_data <- extract_read_data_hdf5r(file_path)
    sampling_rate <- read_data$sampling_rate

    # Z-normalize the data
    norm_data <- z_normalize(read_data$raw_data)

    # Recitfy the data (flip anything that is below zero)
    rectified_data <- rectify(norm_data)

    # Windsorize the data (clip anything above a threshold)
    truncated_data <- truncate_spikes(rectified_data,
                                      spike_threshold=POLY_T_CNDA_SPIKE_THRESHOLD)

    # smoothen the data
    smoothed_data_1 <- left_to_right_sliding_window_cdna_polyt('mean', truncated_data, POLY_T_CDNA_MOVING_WINDOW_SIZE, 1)
    smoothed_data_2 <- right_to_left_sliding_window_cdna_polyt('mean', truncated_data, POLY_T_CDNA_MOVING_WINDOW_SIZE, 1)
    smoothed_data_3 <- pmin(smoothed_data_1, smoothed_data_2)
    smoothed_data <- smoothed_data_3

    # find intersections with the threshold
    intersections <- smoothed_data < POLY_T_CNDA_THRESHOLD
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

    # find ploy(T) tail
    cdna_poly_t_read_type <- ''

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
            ggplot2::geom_hline(yintercept=POLY_T_CNDA_THRESHOLD, color = "black") +
            ggplot2::scale_x_continuous(limits = c(1,
                                                   ceiling(length(smoothed_data_1)/1)))

        if (exists('pri_poly_t_start')) {
            p <- p+ggplot2::geom_line(ggplot2::aes(y = c(rep(NA, times=pri_poly_t_start-1),
                                                         truncated_data[pri_poly_t_start:pri_poly_t_end],
                                                         rep(NA, times=(read_length-pri_poly_t_end)))), color='red')+
                ggplot2::geom_line(ggplot2::aes(y = smoothed_data_3), color='orange')
        }

        if (exists('sec1_poly_t_start')) {
            p <- p+ggplot2::geom_line(ggplot2::aes(y = c(rep(NA, times=sec1_poly_t_start-1),
                                                         truncated_data[sec1_poly_t_start:sec1_poly_t_end],
                                                         rep(NA, times=(read_length-sec1_poly_t_end)))), color='red')+
                ggplot2::geom_line(ggplot2::aes(y = smoothed_data_3), color='orange')
        }

        if (exists('sec2_poly_t_start')) {
            p <- p+ggplot2::geom_line(ggplot2::aes(y = c(rep(NA, times=sec2_poly_t_start-1),
                                                         truncated_data[sec2_poly_t_start:sec2_poly_t_end],
                                                         rep(NA, times=(read_length-sec2_poly_t_end)))), color='red')+
                ggplot2::geom_line(ggplot2::aes(y = smoothed_data_3), color='orange')
        }

        if (show_plots){
            (p)
        }
        p
        if (save_plots){
            filename <- basename(file_path)
            filename <- paste(filename,'.png', sep='')
            dir.create(file.path(save_dir, 'plots', fsep = .Platform$file.sep))
            save_path <- file.path(save_dir, 'plots', filename, fsep = .Platform$file.sep)
            ggplot2::ggsave(save_path, width = 300, height = 70, units = 'mm')
        }

    }
}

