#' Find Poly(A) tail in a single RNA read.
#'
#' @param file_path Path of the FAST5 file
#'
#' @return A list of Fast5 file data
#' @export
#'
#' @examples
#' find_rna_polya_tail_per_read('path/to/fast5/file')
find_rna_polya_tail_per_read <- function(file_path,
                                         save_plots=FALSE,
                                         show_plots=FALSE,
                                         save_dir='~',
                                         plotting_library='rbokeh',
                                         plot_debug=FALSE) {
    print(file_path)
    # Empirical parameters
    POLY_A_RNA_THRESHOLD <- 0.3
    POLY_A_RNA_SPIKE_THRESHOLD <- 3.0
    POLY_A_RNA_MOVING_WINDOW_SIZE <- 400
    POLY_A_RNA_ADAPTOR_TROUGH_SIZE_THRESHOLD <- 2200

    # Read the FAST5 data
    read_data <- extract_read_data_hdf5r(file_path)
    sampling_rate <- read_data$sampling_rate

    # Z-normalize the data
    norm_data <- z_normalize(read_data$raw_data)

    # Recitfy the data (flip anything that is below zero)
    rectified_data <- norm_data

    # Windsorize the data (clip anything above a threshold)
    truncated_data <- truncate_spikes(rectified_data,
                                      spike_threshold=POLY_A_RNA_SPIKE_THRESHOLD)

    # Smoothen the data
    smoothed_data_1 <- left_to_right_sliding_window_cdna_polya('mean', truncated_data, POLY_A_RNA_MOVING_WINDOW_SIZE, 1)
    smoothed_data_2 <- right_to_left_sliding_window_cdna_polya('mean', truncated_data, POLY_A_RNA_MOVING_WINDOW_SIZE, 1)
    smoothed_data_3 <- pmax(smoothed_data_1, smoothed_data_2)
    smoothed_data <- smoothed_data_3

    # Find intersections with the threshold
    intersections <- smoothed_data > POLY_A_RNA_THRESHOLD
    rle_intersections <- rle(intersections)

    # Smoothen RLE intersecitons to remove the fucking W pattern formed by the two overlapping windows
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
                                                                save_plots, show_plots)
    len_rle <- length(rle_values)
    read_length <- length(read_data$raw_data)

    # Find moves in post-poly(A) region
    moves_in_non_poly_a_region <-  NA
    non_poly_a_seq_start <- NA
    non_poly_a_seq_end <- NA
    poly_a_fastq <- NA
    poly_a_length_in_nucleotides_1 <- NA
    poly_a_length_in_nucleotides_2 <- NA

    if (!is.na(precise_polya_boundries$start)) {
        poly_a_fastq <- extract_fastq_in_interval(read_data$event_data,
                                                  precise_polya_boundries$start,
                                                  precise_polya_boundries$end)
        non_poly_a_seq_start <- precise_polya_boundries$end + 1
        non_poly_a_seq_end <- read_length

        moves_in_non_poly_a_region <- extract_moves_in_interval(read_data$event_data,
                                                                non_poly_a_seq_start,
                                                                non_poly_a_seq_end)

        poly_a_length_in_nucleotides_1 <- ((precise_polya_boundries$end - precise_polya_boundries$start) * moves_in_non_poly_a_region) / (non_poly_a_seq_end - non_poly_a_seq_start)

        non_polya_length <- find_non_polya_region_length_without_homopolymers(non_poly_a_seq_start,
                                                                              non_poly_a_seq_end,
                                                                              read_data$moves_sample_wise_vector,
                                                                              homopolymer_threshold = 50)

        poly_a_length_in_nucleotides_2 <- ((precise_polya_boundries$end - precise_polya_boundries$start) * moves_in_non_poly_a_region) / non_polya_length
    }

    if (show_plots || save_plots) {
        filename <- basename(file_path)

        df = data.frame(x=c(1:length(read_data$raw_data)),
                        truncated_data=truncated_data,
                        smoothed_data_1=smoothed_data_1,
                        smoothed_data_2=smoothed_data_2,
                        smoothed_data_3=smoothed_data_3,
                        moves=read_data$moves_sample_wise_vector,
                        mean_data=precise_polya_boundries$mean_data,
                        slope=precise_polya_boundries$slope)
        # add a poly(A) tail to the dataframe if it exists
        if (!is.na(precise_polya_boundries$start)) {
            df['poly_a_tail'] <- c(rep(NA, times=precise_polya_boundries$start-1),
                                   truncated_data[precise_polya_boundries$start:precise_polya_boundries$end],
                                   rep(NA, times=(read_length-precise_polya_boundries$end)))
        }

        if (plotting_library == 'ggplot2') {
            p <- ggplot2::ggplot(data=df, ggplot2::aes(x = x)) +
                ggplot2::geom_line(ggplot2::aes(y = truncated_data), color='blue')
            if (plot_debug) {
                p <- p + ggplot2::geom_line(ggplot2::aes(y = moves)) +
                ggplot2::geom_hline(yintercept=POLY_A_RNA_THRESHOLD, color = "black") +
                ggplot2::geom_hline(yintercept=-POLY_A_RNA_THRESHOLD, color = "black") +
                ggplot2::scale_x_continuous(limits = c(1, ceiling(length(smoothed_data_1)/1)))
            }
            if (!is.na(precise_polya_boundries$start)) {
                p <- p + ggplot2::geom_line(ggplot2::aes(y = c(rep(NA, times=precise_polya_boundries$start-1),
                                                               truncated_data[precise_polya_boundries$start:precise_polya_boundries$end],
                                                               rep(NA, times=(read_length-precise_polya_boundries$end)))), color='red') +
                ggplot2::scale_colour_manual(values =c('red'='red', 'blue'='blue'), labels = c('Poly(A) tail','Normalized windsorized data'))
            }
            if (plot_debug){
                p <- p + ggplot2::geom_line(ggplot2::aes(y = smoothed_data_3), color='orange')
            }
            p <- p + ggplot2::ggtitle(filename) +
                ggplot2::xlab('Sample index') +
                ggplot2::ylab('Z-normalized data values')
        }

        if (plotting_library == 'rbokeh') {
            p <- rbokeh::figure(data=df, width=1000, height=600, title=filename)
            p <- rbokeh::ly_lines(p, truncated_data, color='blue', legend = "Normalized windsorized data")
            if (plot_debug) {
                p <- rbokeh::ly_lines(p, slope, color='red', width=3, legend = "Slope")
                p <- rbokeh::ly_lines(p, mean_data, color='green', width=3, legend = "Mean of potential poly(A) region")
                p <- rbokeh::ly_lines(p, moves, color='black', width=1, legend = "Moves")
                p <- rbokeh::ly_lines(p, smoothed_data_3, color='orange', width=3, legend = "Smoothed data")
                p <- rbokeh::ly_abline(p, h=POLY_A_RNA_THRESHOLD, color = 'black', type = 2, width=2, legend = "Poly(A) threshold")
                p <- rbokeh::ly_abline(p, h=-POLY_A_RNA_THRESHOLD, color = 'black',type = 2,  width=2, legend = "Poly(A) threshold")
                if (!is.na(precise_polya_boundries$start)) {
                    p <- rbokeh::ly_abline(p, v=precise_polya_boundries$start, color = 'orange', type = 2, width=2, legend = "Poly(A) start")
                    p <- rbokeh::ly_abline(p, v=precise_polya_boundries$end, color = 'red', type = 2, width=2, legend = "Poly(A) end")
                }
            } else {
                p <- rbokeh::ly_lines(p, poly_a_tail, color = 'red', legend = "Poly(A) tail")
            }
            p <- rbokeh::x_axis(p, label='Sample index')
            p <- rbokeh::y_axis(p, label='Z-normalized data values')
            p <- rbokeh::tool_pan(p, dimensions = "width")
            p <- rbokeh::tool_wheel_zoom(p, dimensions = "width")
        }

        if (save_plots) {
            dir.create(file.path(save_dir, 'plots', fsep = .Platform$file.sep))
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

    data <- list(read_id=read_data$read_id,
                 poly_a_start=precise_polya_boundries$start,
                 poly_a_end=precise_polya_boundries$end,
                 poly_a_fastq=poly_a_fastq,
                 poly_a_length_in_nucleotides_1=poly_a_length_in_nucleotides_1,
                 poly_a_length_in_nucleotides_2=poly_a_length_in_nucleotides_2,

                 non_poly_a_seq_start = non_poly_a_seq_start,
                 non_poly_a_seq_end = non_poly_a_seq_end,
                 moves_in_non_poly_a_region = moves_in_non_poly_a_region,

                 sampling_rate=sampling_rate,
                 file_path=file_path)
    return(data)
}

