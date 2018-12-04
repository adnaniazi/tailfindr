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
                                          poly_a_adaptor="GAAGATAGAGCGACAGGCAAGT",
                                          save_plots=FALSE,
                                          show_plots=FALSE,
                                          save_dir='~'){

    # Empirical parameters
    POLY_T_CNDA_THRESHOLD <- 0.20
    POLY_T_CNDA_SPIKE_THRESHOLD <- 2.0
    POLY_T_CDNA_MOVING_WINDOW_SIZE <- 120

    # Read the FAST5 data
    read_data <- extract_read_data_hdf5r(file_path)
    sampling_rate <- read_data$sampling_rate

    # Empirical parameters
    # Allow for 40 nt gap only
    POLY_T_CDNA_SEC_POLY_T_MAX_GAP <- read_data$samples_per_nt * 40

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

    cdna_poly_t_read_type <- ''
    non_poly_t_seq_start <- NA
    non_poly_t_seq_end <- NA
    moves_in_non_poly_t_region <- NA

    # Do  a second pass with smaller windows to capture shorter tail that the wider window migth have missed
    if (length(rle_values) <= 2) {
        # smoothen the data
        POLY_T_CDNA_MOVING_WINDOW_SIZE <- POLY_T_CDNA_MOVING_WINDOW_SIZE/4
        POLY_T_CNDA_THRESHOLD <- POLY_T_CNDA_THRESHOLD
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
    }


    if (rle_values[1]){
        # ignore the bullshit poly(T) looking stuff at the beginning
        rle_start <- 2
        cdna_poly_t_read_type_begin <- '1'
    } else {
        rle_start <- 1
        cdna_poly_t_read_type_begin <- ''
    }

    # find ploy(T) tail
    # CASE 1: bullshit --> adaptor --> tail --> some_sequence
    # file: 0.fast5
    #$lengths
    #[1]  719  399  400 8566
    #$values
    #[1]  TRUE FALSE  TRUE FALSE

    if (len_rle >= 3 && !rle_values[rle_start] && rle_values[rle_start+1]) {
        pri_poly_t_start <- rle_indices[rle_start]
        pri_poly_t_end <- rle_indices[rle_start+1]
        pri_poly_t_fastq <- extract_fastq_in_interval(read_data$event_data,
                                                      pri_poly_t_start,
                                                      pri_poly_t_end)
        cdna_poly_t_read_type <- paste(cdna_poly_t_read_type_begin, 'adaptor10', sep = '')
        non_poly_t_seq_start <- pri_poly_t_end + 1
        non_poly_t_seq_end <- rle_indices[rle_start+2]

        # find the first secondary tail
        if (len_rle >= 5 && !rle_values[rle_start+2] && rle_values[rle_start+3] &&
            rle_lengths[rle_start+2] < POLY_T_CDNA_SEC_POLY_T_MAX_GAP) {
            gap1_start <- pri_poly_t_end + 1
            gap1_end <- rle_indices[rle_start+2] - 1
            sec1_poly_t_start <- rle_indices[rle_start+2]
            sec1_poly_t_end <- rle_indices[rle_start+3]
            sec1_poly_t_fastq <- extract_fastq_in_interval(read_data$event_data,
                                                           sec1_poly_t_start,
                                                           sec1_poly_t_end)
            gap1_fastq <- extract_fastq_in_interval(read_data$event_data,
                                                    gap1_start,
                                                    gap1_end)
            cdna_poly_t_read_type <- paste(cdna_poly_t_read_type_begin, '10', sep = '')
            non_poly_t_seq_start <- sec1_poly_t_end + 1
            non_poly_t_seq_end <-  rle_indices[rle_start+4]

            # find the second secondary tail
            if (len_rle >= 7 && !rle_values[rle_start] && rle_values[rle_start+1] &&
                rle_lengths[rle_start+4] < POLY_T_CDNA_SEC_POLY_T_MAX_GAP) {
                gap2_start <- sec1_poly_t_end + 1
                gap2_end <- rle_indices[rle_start+4] - 1
                sec2_poly_t_start <- rle_indices[rle_start+4]
                sec2_poly_t_end <- rle_indices[rle_start+5]
                sec2_poly_t_fastq <- extract_fastq_in_interval(read_data$event_data,
                                                               sec2_poly_t_start,
                                                               sec2_poly_t_end)
                gap2_fastq <- extract_fastq_in_interval(read_data$event_data,
                                                        gap2_start,
                                                        gap2_end)
                cdna_poly_t_read_type <- paste(cdna_poly_t_read_type_begin, '10', sep = '')
                non_poly_t_seq_start <- sec2_poly_t_end + 1
                non_poly_t_seq_end <-  rle_indices[rle_start+6]
            }
        }
    }

    # find the adaptor attached to the beginning of primary poly(T) tail by
    # aligning anything that comes before the primary poly(T) tail to the ONT-provided adaptor sequences
    # adaptor seq: ACTTGCCTGTCGCTCTATCTTC | 22
    if (exists('pri_poly_t_start')){
        if (pri_poly_t_start == 0){
            tail_adaptor <- paste('Tail adaptor absent; aln score: NA; adaptor seq: NA')
            has_valid_poly_t_tail <- FALSE
        } else {
            ta_hvptt <- align_cdna_polyt_adaptor(read_data$event_data, pri_poly_t_start, poly_a_adaptor)
            tail_adaptor_seq <- ta_hvptt$tail_adaptor_seq
            tail_adaptor_aln_score <- ta_hvptt$tail_adaptor_aln_score
            has_valid_poly_t_tail <- ta_hvptt$has_valid_poly_t_tail
        }

    } else {
        tail_adaptor <- NA
        has_valid_poly_t_tail <- FALSE
    }

    # Calculate the length of non-poly-T region
    # This information can be used to find the dwell time per basepair
    if (!is.na(non_poly_t_seq_start)) {
        if (is.na(non_poly_t_seq_end)) {
            non_poly_t_seq_end <- length(norm_data)
        }
        moves_in_non_poly_t_region <- extract_moves_in_interval(read_data$event_data, non_poly_t_seq_start, non_poly_t_seq_end)
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
            print(p)
        }
        if (save_plots){
            filename <- basename(file_path)
            filename <- paste(filename,'.png', sep='')
            dir.create(file.path(save_dir, 'plots', fsep = .Platform$file.sep))
            save_path <- file.path(save_dir, 'plots', filename, fsep = .Platform$file.sep)
            ggplot2::ggsave(save_path, plot = p, width = 300, height = 70, units = 'mm')
        }
    }

    if (!exists('sec2_poly_t_start')){
        sec2_poly_t_start <- NA
        sec2_poly_t_end <- NA
        sec2_poly_t_fastq <- NA
        gap2_start <- NA
        gap2_end <- NA
        gap2_fastq <- NA
    }

    if (!exists('sec1_poly_t_start')){
        sec1_poly_t_start <- NA
        sec1_poly_t_end <- NA
        sec1_poly_t_fastq <- NA
        gap1_start <- NA
        gap1_end <- NA
        gap1_fastq <- NA
    }

    if (!exists('pri_poly_t_start')){
        pri_poly_t_start <- NA
        pri_poly_t_end <- NA
        cdna_poly_t_read_type <- 'Tail not found'
        has_valid_poly_t_tail <- FALSE
        pri_poly_t_fastq <- NA
    }

    data <- list(read_id = read_data$read_id,

                 pri_poly_t_start = pri_poly_t_start,
                 pri_poly_t_end = pri_poly_t_end,
                 pri_poly_t_fastq = pri_poly_t_fastq,

                 gap1_start = gap1_start,
                 gap1_end = gap1_end,
                 gap1_fastq = gap1_fastq,

                 sec1_poly_t_start = sec1_poly_t_start,
                 sec1_poly_t_end = sec1_poly_t_end,
                 sec1_poly_t_fastq = sec1_poly_t_fastq,

                 gap2_start = gap2_start,
                 gap2_end = gap2_end,
                 gap2_fastq = gap2_fastq,

                 sec2_poly_t_start = sec2_poly_t_start,
                 sec2_poly_t_end = sec2_poly_t_end,
                 sec2_poly_t_fastq = sec2_poly_t_fastq,

                 non_poly_t_seq_start = non_poly_t_seq_start,
                 non_poly_t_seq_end = non_poly_t_seq_end,
                 moves_in_non_poly_t_region = moves_in_non_poly_t_region,

                 sampling_rate = sampling_rate,
                 cdna_poly_t_read_type = cdna_poly_t_read_type,
                 tail_adaptor_seq = tail_adaptor_seq,
                 tail_adaptor_aln_score = tail_adaptor_aln_score,
                 has_valid_poly_t_tail = has_valid_poly_t_tail,
                 samples_per_nt = read_data$samples_per_nt,
                 file_path = file_path)
    return(data)
}

