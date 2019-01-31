#' Find Poly(A) tail in a single RNA read.
#'
#' @param file_path Path of the FAST5 file
#' @param read_id_fast5_file
#' @param multifast5
#' @param basecalled_with
#' @param save_plots
#' @param show_plots
#' @param save_dir
#' @param plotting_library
#' @param plot_debug
#'
#' @return A list of Fast5 file data
#' @export
#'
#' @examples
#' find_rna_polya_tail_per_read('path/to/fast5/file')
find_rna_polya_tail_per_read <- function(file_path = NA,
                                         read_id_fast5_file = NA,
                                         multifast5,
                                         basecalled_with,
                                         save_plots = FALSE,
                                         show_plots = FALSE,
                                         save_dir = '~',
                                         plotting_library = 'rbokeh',
                                         plot_debug = FALSE) {
    # Empirical parameters
    POLY_A_RNA_THRESHOLD <- 0.3
    POLY_A_RNA_SPIKE_THRESHOLD <- 3.0
    POLY_A_RNA_MOVING_WINDOW_SIZE <- 400
    POLY_A_RNA_ADAPTOR_TROUGH_SIZE_THRESHOLD <- 2200

    # Read the FAST5 data
    read_data <- extract_read_data(file_path,
                                   read_id_fast5_file,
                                   plot_debug,
                                   basecalled_with,
                                   multifast5,
                                   model)
    # first read the data and find the tailtype
    if (multifast5) {
        file_path <- read_id_fast5_file$fast5_file
    }
    read_id <- read_data$read_id

    # Z-normalize the data
    norm_data <- z_normalize(read_data$raw_data)

    # Recitfy the data (flip anything that is below zero)
    rectified_data <- norm_data

    # Windsorize the data (clip anything above a threshold)
    truncated_data <- truncate_spikes(rectified_data,
                                      spike_threshold=POLY_A_RNA_SPIKE_THRESHOLD)

    # Smoothen the data
    smoothed_data_1 <- left_to_right_sliding_window_cdna_polya('mean',
                                                               truncated_data,
                                                               POLY_A_RNA_MOVING_WINDOW_SIZE, 1)
    smoothed_data_2 <- right_to_left_sliding_window_cdna_polya('mean',
                                                               truncated_data,
                                                               POLY_A_RNA_MOVING_WINDOW_SIZE, 1)
    smoothed_data_3 <- pmax(smoothed_data_1, smoothed_data_2)
    smoothed_data <- smoothed_data_3

    # Find intersections with the threshold
    intersections <- smoothed_data > POLY_A_RNA_THRESHOLD
    rle_intersections <- rle(intersections)

    # Smoothen RLE intersecitons to remove the annoying W pattern formed by the
    # merging the two smoothened signals that have a relative shift with each
    # other.
    rle_intersections <- smoothen_rle_intersections_rna(rle_intersections)
    rle_lengths <- rle_intersections$lengths
    rle_values <- rle_intersections$values
    rle_indices <- cumsum(rle_lengths)

    # Find crude poly(A) boundries
    crude_polya_boundries <- find_rna_polya_crude_start_end(rle_lengths,
                                                            rle_values,
                                                            rle_indices,
                                                            POLY_A_RNA_ADAPTOR_TROUGH_SIZE_THRESHOLD)

    # Refine the curde poly(A) boundries
    precise_polya_boundries <- find_rna_polya_precise_start_end(truncated_data,
                                                                POLY_A_RNA_THRESHOLD,
                                                                crude_polya_boundries,
                                                                save_plots,
                                                                show_plots)
    len_rle <- length(rle_values)
    read_length <- length(read_data$raw_data)

    # Find moves in the poly(A) region
    poly_a_fastq <- NA
    tail_length <- NA
    tail_start <- NA
    tail_end <- NA
    if (!is.na(precise_polya_boundries$start) & !is.na(precise_polya_boundries$end)) {
        poly_a_fastq <- extract_fastq_in_interval(read_data$event_data,
                                                  precise_polya_boundries$start,
                                                  precise_polya_boundries$end)
        tail_start <- precise_polya_boundries$start
        tail_end <- precise_polya_boundries$end
        tail_length <- (tail_end - tail_start) / read_data$samples_per_nt
    }

    if (show_plots | save_plots) {
        if (multifast5){
            filename <- read_id
        } else {
            filename <- basename(file_path)
        }
        plot_title <- paste('Poly(A) tail  |  ',
                            'Tail length[nt]: ', round(tail_length, 2), '  |  ',
                            'Tail start: ', tail_start, '  |  ',
                            'Tail end: ', tail_end, '  |  ',
                            'Tail duration[Sa]: ', tail_end-tail_start, '  |  ',
                            'Samples per nt: ', round(read_data$samples_per_nt, 2),
                            sep='')

        df = data.frame(x=c(1:length(read_data$raw_data)),
                        truncated_data=truncated_data,
                        smoothed_data_3=smoothed_data_3,
                        moves=read_data$moves_sample_wise_vector + 1,
                        mean_data=precise_polya_boundries$mean_data,
                        slope=precise_polya_boundries$slope)

        # add a poly(A) tail to the dataframe if it exists
        if (!is.na(tail_start) & !is.na(tail_end)) {
            df['poly_a_tail'] <- c(rep(NA, times=tail_start-1),
                                   truncated_data[tail_start:tail_end],
                                   rep(NA, times=(read_length-tail_end)))
        }

        if (plotting_library == 'ggplot2') {
            p <- ggplot2::ggplot(data=df, ggplot2::aes(x = x))
            if (plot_debug) {
                p <- p + ggplot2::geom_line(ggplot2::aes(y = truncated_data), color='#4040a1') +
                ggplot2::geom_line(ggplot2::aes(y = moves), color = '#b2b2b2') +
                ggplot2::geom_hline(yintercept = POLY_A_RNA_THRESHOLD, color = "#00B0DF") +
                ggplot2::geom_line(ggplot2::aes(y = slope, color = "#BF1268")) +
                ggplot2::geom_hline(yintercept = -POLY_A_RNA_THRESHOLD, color = "#FF94CD")
            } else {
                p <- p + ggplot2::geom_line(ggplot2::aes(y = truncated_data), color='#8E8E8E')
            }

            if (!is.na(tail_start) & !is.na(tail_end)) {
                p <- p + ggplot2::geom_line(ggplot2::aes(y = c(rep(NA, times = tail_start-1),
                                                               truncated_data[tail_start:tail_end],
                                                               rep(NA, times = (read_length-tail_end)))), color='#BF1268')
            }
            if (plot_debug){
                p <- p + ggplot2::geom_line(ggplot2::aes(y = smoothed_data_3), color='#060B54') +
                ggplot2::geom_line(ggplot2::aes(y = mean_data, color = "#F79A14"))
            }
            p <- p + ggplot2::ggtitle(plot_title) +
                ggplot2::xlab('Sample index') +
                ggplot2::ylab('z-normalized data')
        }

        if (plotting_library == 'rbokeh') {
            p <- rbokeh::figure(data=df, width=1000, height=600, title=plot_title)
            p <- rbokeh::ly_lines(p, truncated_data, color='#4040a1', legend = "Normalized windsorized data")
            if (plot_debug) {
                p <- rbokeh::ly_lines(p, slope, color='#BF1268', width=3, legend = "Slope")
                p <- rbokeh::ly_lines(p, mean_data, color='#F79A14', width = 3, legend = "Mean of potential poly(A) region")
                p <- rbokeh::ly_lines(p, moves, color='#b2b2b2', width=1, legend = "Moves")
                p <- rbokeh::ly_lines(p, smoothed_data_3, color='#060B54', width=3, legend = "Smoothed data")
                p <- rbokeh::ly_abline(p, h=POLY_A_RNA_THRESHOLD, color = '#00B0DF', type = 3, width=2, legend = "Slope upper bound")
                p <- rbokeh::ly_abline(p, h=-POLY_A_RNA_THRESHOLD, color = '#FF94CD', type = 3,  width=2, legend = "Slope lower bound")
                if (!is.na(precise_polya_boundries$start)) {
                    p <- rbokeh::ly_abline(p, v=precise_polya_boundries$start, color = 'orange', type = 2, width = 2, legend = "Poly(A) start")
                    p <- rbokeh::ly_abline(p, v=precise_polya_boundries$end, color = 'red', type = 2, width = 2, legend = "Poly(A) end")
                }
            } else {
                p <- rbokeh::ly_lines(p, poly_a_tail, color = 'red', legend = "Poly(A) tail")
            }
            p <- rbokeh::x_axis(p, label='Sample index')
            p <- rbokeh::y_axis(p, label='z-normalized data')
            p <- rbokeh::tool_pan(p, dimensions = "width")
            p <- rbokeh::tool_wheel_zoom(p, dimensions = "width")
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

        if (show_plots) {
            print(p)
        }
    }

    data <- list(read_id = read_id,
                 tail_start = precise_polya_boundries$start,
                 tail_end = precise_polya_boundries$end,
                 samples_per_nt = read_data$samples_per_nt,
                 tail_length = tail_length,
                 polya_fastq = poly_a_fastq,
                 file_path = file_path)
    return(data)
}

