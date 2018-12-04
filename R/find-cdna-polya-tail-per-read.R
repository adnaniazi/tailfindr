#' Find Poly(A) tail in a single cDNA read. The read must be a poly(A) read.
#'
#' @param file_path Path of the FAST5 file
#'
#' @return A list of Fast5 file data
#' @export
#'
#' @examples
#' find_cdna_polya_tail_per_read('path/to/fast5/file')
find_cdna_polya_tail_per_read <- function(file_path,
                                          poly_a_adaptor="GAAGATAGAGCGACAGGCAAGT",
                                          save_plots=FALSE,
                                          show_plots=FALSE,
                                          save_dir='~'){

    # read the FAST5 data
    read_data <- extract_read_data_hdf5r(file_path)
    sampling_rate <- read_data$sampling_rate

    # Empirical parameters
    POLY_A_CNDA_THRESHOLD <- 0.31
    POLY_A_CNDA_SPIKE_THRESHOLD <- 2.0
    POLY_A_CDNA_MOVING_WINDOW_SIZE <- 120
    POLY_A_CDNA_SEC_POLY_A_MAX_GAP <- 1200
    POLY_A_CDNA_SEC_POLY_A_MIN_SIZE <- read_data$samples_per_nt * 10

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
    non_poly_a_seq_start <- NA
    non_poly_a_seq_end <- NA
    moves_in_non_poly_a_region <- NA

    # case X
    # file: 1.fast5
    # Run Length Encoding
    # lengths: int [1:2] 946 14784
    # values : logi [1:2] TRUE FALSE
    if (len_rle <= 2){
        # Do  a second pass with smaller windows to capture shorter tail that the wider window migth have missed
        cdna_poly_a_read_type <- 'second-round'
        POLY_A_CDNA_MOVING_WINDOW_SIZE <- POLY_A_CDNA_MOVING_WINDOW_SIZE/4
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
    }

    # case X
    # file: 5.fast5, 0.fast5
    # from right to left
    # bullshit(FALSE) --> pri_polyA(TRUE) --> bullshit(FALSE) of some length> POLY_A_CDNA_SEC_POLY_A_MAX_GAP
    # $lengths: 1901 9877  734  360
    # $values:TRUE FALSE  TRUE FALSE
    if (!rle_values[len_rle] && rle_values[(len_rle-1)] &&
             !rle_values[(len_rle-2)]){
        pri_poly_a_start <- rle_indices[(len_rle-2)]
        pri_poly_a_end <- rle_indices[(len_rle-1)]
        pri_poly_a_fastq <- extract_fastq_in_interval(read_data$event_data, pri_poly_a_start, pri_poly_a_end)
        cdna_poly_a_read_type <- '..010'
        if (length(rle_indices)==3){
            non_poly_a_seq_start <- 0
        } else {
            non_poly_a_seq_start <- rle_indices[(len_rle-3)]
        }
        non_poly_a_seq_end <- pri_poly_a_start - 1

        # find the first secondary tail
        if (len_rle > 4) {
            # if we have a small gap then we a have the first secondary poly-a tail
            # file: 0.fast5
            if ((rle_lengths[(len_rle-2)] < POLY_A_CDNA_SEC_POLY_A_MAX_GAP) &&
                (rle_lengths[(len_rle-3)] > POLY_A_CDNA_SEC_POLY_A_MIN_SIZE) &&
                rle_values[(len_rle-3)] && !rle_values[(len_rle-4)]){
                sec1_poly_a_start <- rle_indices[(len_rle-4)]
                sec1_poly_a_end <- rle_indices[(len_rle-3)]
                sec1_poly_a_fastq <- extract_fastq_in_interval(read_data$event_data, sec1_poly_a_start, sec1_poly_a_end)
                gap1_start <- sec1_poly_a_end + 1
                gap1_end <- pri_poly_a_start - 1
                gap1_fastq <- extract_fastq_in_interval(read_data$event_data, gap1_start, gap1_end)
                cdna_poly_a_read_type <- '..01010'
                if (length(rle_indices)==5){
                    non_poly_a_seq_start <- 0
                } else {
                    non_poly_a_seq_start <- rle_indices[(len_rle-5)]
                }
                non_poly_a_seq_end <- sec1_poly_a_start - 1

                # if first secondary tails is found, then find the second secondary tail
                if (len_rle > 6) {
                    # if we have a small gap then we a have the first secondary poly-a tail
                    # file: 0.fast5
                    if ((rle_lengths[(len_rle-4)] < POLY_A_CDNA_SEC_POLY_A_MAX_GAP)&&
                        (rle_lengths[(len_rle-5)] > POLY_A_CDNA_SEC_POLY_A_MIN_SIZE) &&
                        rle_values[(len_rle-5)] && !rle_values[(len_rle-6)]){
                        sec2_poly_a_start <- rle_indices[(len_rle-6)]
                        sec2_poly_a_end <- rle_indices[(len_rle-5)]
                        sec2_poly_a_fastq <- extract_fastq_in_interval(read_data$event_data, sec2_poly_a_start, sec2_poly_a_end)
                        gap2_start <- sec2_poly_a_end + 1
                        gap2_end <- sec1_poly_a_start - 1
                        gap2_fastq <- extract_fastq_in_interval(read_data$event_data, gap2_start, gap2_end)
                        cdna_poly_a_read_type <- '..0101010'
                        if (length(rle_indices)==7){
                            non_poly_a_seq_start <- 0
                        } else {
                            non_poly_a_seq_start <- rle_indices[(len_rle-7)]
                        }
                        non_poly_a_seq_end <- sec2_poly_a_start - 1
                    }
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
        pri_poly_a_fastq <- extract_fastq_in_interval(read_data$event_data, pri_poly_a_start, pri_poly_a_end)
        cdna_poly_a_read_type <- 'adaptor01'
        non_poly_a_seq_start <- rle_indices[(len_rle-2)]
        non_poly_a_seq_end <- pri_poly_a_start - 1
    }

    # find the adaptor attached to the end of primary poly(A) tail by
    # aligning anything that comes after the primary poly(A) tail to the ONT-provided adaptor sequences
    # adaptor seq: GAAGATAGAGCGACAGGCAAGT | 22
    if (exists('pri_poly_a_start')){
        if (pri_poly_a_end == read_length){
            tail_adaptor <- paste('Tail adaptor absent; aln score: NA; adaptor seq: NA')
            has_valid_poly_a_tail <- FALSE
        } else {
            ta_hvpat <- align_cdna_polya_adaptor(read_data$event_data, pri_poly_a_end, poly_a_adaptor)
            tail_adaptor_seq <- ta_hvpat$tail_adaptor_seq
            tail_adaptor_aln_score <- ta_hvpat$tail_adaptor_aln_score
            has_valid_poly_a_tail <- ta_hvpat$has_valid_poly_a_tail
        }

    } else {
        tail_adaptor <- NA
        has_valid_poly_a_tail <- FALSE
    }

    # Calculate the length of non-poly-A region
    # This information can be used to find the dwell time per basepair
    if (!is.na(non_poly_a_seq_start)) {
        moves_in_non_poly_a_region <- extract_moves_in_interval(read_data$event_data, non_poly_a_seq_start, non_poly_a_seq_end)
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

        if (show_plots){print(p)}
        if (save_plots){
            filename <- basename(file_path)
            filename <- paste(filename,'.png', sep='')
            dir.create(file.path(save_dir, 'plots', fsep = .Platform$file.sep))
            save_path <- file.path(save_dir, 'plots', filename, fsep = .Platform$file.sep)
            ggplot2::ggsave(save_path, plot = p, width = 300, height = 70, units = 'mm')
        }

    }

    if (!exists('sec2_poly_a_start')){
        sec2_poly_a_start <- NA
        sec2_poly_a_end <- NA
        sec2_poly_a_fastq <- NA
        gap2_start <- NA
        gap2_end <- NA
        gap2_fastq <- NA
    }

    if (!exists('sec1_poly_a_start')){
        sec1_poly_a_start <- NA
        sec1_poly_a_end <- NA
        sec1_poly_a_fastq <- NA
        gap1_start <- NA
        gap1_end <- NA
        gap1_fastq <- NA
    }

    if (!exists('pri_poly_a_start')){
        pri_poly_a_start <- NA
        pri_poly_a_end <- NA
        cdna_poly_a_read_type <- 'Tail not found'
        has_valid_poly_a_tail <- FALSE
        pri_poly_a_fastq <- NA
    }

    data <- list(read_id=read_data$read_id,

                 pri_poly_a_start=pri_poly_a_start,
                 pri_poly_a_end=pri_poly_a_end,
                 pri_poly_a_fastq=pri_poly_a_fastq,

                 gap1_start=gap1_start,
                 gap1_end=gap1_end,
                 gap1_fastq=gap1_fastq,

                 sec1_poly_a_start=sec1_poly_a_start,
                 sec1_poly_a_end=sec1_poly_a_end,
                 sec1_poly_a_fastq=sec1_poly_a_fastq,

                 gap2_start=gap2_start,
                 gap2_end=gap2_end,
                 gap2_fastq=gap2_fastq,

                 sec2_poly_a_start=sec2_poly_a_start,
                 sec2_poly_a_end=sec2_poly_a_end,
                 sec2_poly_a_fastq=sec2_poly_a_fastq,

                 non_poly_a_seq_start = non_poly_a_seq_start,
                 non_poly_a_seq_end = non_poly_a_seq_end,
                 moves_in_non_poly_a_region = moves_in_non_poly_a_region,

                 sampling_rate=sampling_rate,
                 cdna_poly_a_read_type=cdna_poly_a_read_type,
                 tail_adaptor_seq = tail_adaptor_seq,
                 tail_adaptor_aln_score = tail_adaptor_aln_score,
                 has_valid_poly_a_tail=has_valid_poly_a_tail,
                 samples_per_nt = read_data$samples_per_nt,
                 file_path=file_path)
    return(data)
}

