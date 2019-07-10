#' Title
#'
#' @param rca_data
#' @param read_data
#' @param file_path
#'
#' @return
#' @export
#'
#' @importFrom magrittr "%>%"
#'
#' @examples
rca_find_tails_per_read <- function(rca_data, read_data, file_path) {

    if (is.null(rca_data)) {
        return(data.frame(
            start = NA,
            end = NA,
            fastq_segment = NA,
            read_type = 'invalid',
            cluster = NA,
            tail_start = NA,
            tail_end = NA,
            read_id = read_data$read_id,
            cluster_n = NA,
            file_path = file_path,
            tail_length = NA,
            samples_per_nt = NA
        ))
    }

    # Define thresholds
    SPIKE_THRESHOLD <- 2.0
    MOVING_WINDOW_SIZE <- 30
    MAX_GAP_BETWEEN_TAILS <- 90 * read_data$samples_per_nt
    SLOPE_THRESHOLD <- 0.20

    # Get raw data and event data
    raw_data <- read_data$raw_data
    event_data <- read_data$event_data
    samples_per_nt <- read_data$samples_per_nt
    raw_data_len <- length(raw_data)

    # Normalize data
    norm_data <- z_normalize(read_data$raw_data)

    # Rectifiy data
    rectified_data <- rectify(norm_data)

    # Windsorize the data (clip anything above a threshold)
    truncated_data_org <- truncate_spikes(rectified_data,
                                      spike_threshold = SPIKE_THRESHOLD)

    # filter out consensus data
    read_type <- NULL # for CRAN
    rca_data <- rca_data %>%
        dplyr::mutate(tail_start = NA, tail_end = NA, tail_length = NA) %>%
        dplyr::mutate(read_id = read_data$read_id) %>%
        dplyr::mutate(id = dplyr::row_number())

    rca_data_tmp <- rca_data %>%
        dplyr::filter(read_type == 'polyA' | read_type == 'polyT')

    read_type <- rca_data_tmp$read_type
    fastq_start <- rca_data_tmp$start
    fastq_end <- rca_data_tmp$end
    id <- rca_data_tmp$id

    # Main loop that goes over all the segments and finds poly(A)/poly(T)
    # tail lengths
    for (seg in seq(nrow(rca_data_tmp))) {
        # find bounds to search the tail in
        if (read_type[seg] == 'polyA') {
            bnds <- rca_get_segment_boundaries_in_samples(
                event_data = event_data,
                fastq_start_index = fastq_start[seg],
                fastq_end_index = fastq_end[seg],
                read_type = 'polyA',
                raw_data_len = raw_data_len,
                samples_per_nt = samples_per_nt)
        } else if (read_type[seg] == 'polyT') {
            bnds <- rca_get_segment_boundaries_in_samples(
                event_data = event_data,
                fastq_start_index = fastq_start[seg],
                fastq_end_index = fastq_end[seg],
                read_type = 'polyT',
                raw_data_len = raw_data_len,
                samples_per_nt = samples_per_nt)
        }

        # STEP 1: Prepare data with in the bounds found
        start_index <- bnds$sample_start_idx
        end_index <- bnds$sample_end_idx
        truncated_data <- truncated_data_org[start_index:end_index]
        len_data <-  length(truncated_data)

        # STEP2: Apply smoothening from Left to Right and
        # then from Right to left and then create smoothened data frome the two
        smoothed_data_lr <-
            left_to_right_sliding_window_dna('mean',
                                             truncated_data,
                                             MOVING_WINDOW_SIZE, 1)
        smoothed_data_rl <-
            right_to_left_sliding_window_dna('mean',
                                             truncated_data,
                                             MOVING_WINDOW_SIZE, 1)
        smoothed_data <- pmin(smoothed_data_lr, smoothed_data_rl)

        # STEP X: Reverse poly(A) segment's data to bring it into the
        # same orientation as poly(T) segments
        if (read_type[seg] == 'polyA') {
            # Reverse polyA reads
            truncated_data <- rev(truncated_data)
            smoothed_data <- rev(smoothed_data)
        }

        start <- tail_start <- 1

        # STEP 4: Make slope and mean signals
        window_size <- step_size <- 10
        smd <- rca_make_slope_and_mean_signals(data = truncated_data,
                                               window_size = window_size,
                                               step_size = step_size)
        slope <- smd$slope
        mean_data <- smd$mean_data

        # STEP 5: Refine the left boundary of the tail
        k <- 1
        precise_tail_start <- NA
        while (k < 20) {
            if ((slope[k] < SLOPE_THRESHOLD) & (slope[k] > -SLOPE_THRESHOLD) &
                (mean_data[k] < SLOPE_THRESHOLD + 0.1) & (mean_data[k] > 0)) {
                precise_tail_start <- (k - 1) * window_size + start
                break
            }
            k <- k + 1
        }
        i <- k

        # STEP 6: Refine the right boundary of the tail
        small_glitch_count <- 0
        quit_searching <- FALSE
        tail_end <- NA

        # Main Algorithm:
        # DO while we haven't reached the end of the search window:
        #   - Extend the tail if slope and smoothened data is within the threshold
        #   - if the above criteria fails then:
        #     - check for gap
        #     - there must be at least 60-sample (or more) long tail adjacent to the gap
        #     - the gap should be smaller than 120 nucleotide

        # Define smoothened data threshold first (Different across poly(A)/(T))
        if (read_type[seg] == 'polyA') {
            sm_data_threshold <- 0.6
        } else if  (read_type[seg] == 'polyT') {
            sm_data_threshold <- 0.3
        }

        while (i < length(slope)){
            if ((slope[i] < SLOPE_THRESHOLD) & (slope[i] > -SLOPE_THRESHOLD) &
                (smoothed_data[tail_start+i*window_size] < sm_data_threshold)) {
                tail_end <- i
                i <- i + 1
            } else {
                j <- i
                while (j < length(slope)) {
                    if (j > floor(MAX_GAP_BETWEEN_TAILS/window_size)) {
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
        } # End of while loop

        # tail ends prematurely especially the poly(A) because of the main tail
        # finding logic above. This code extends that tail so that it is correct
        while (i < length(slope)) {
            if (smoothed_data[i*window_size] < sm_data_threshold) {
                tail_end <- i
                i <- i + 1
            } else {
                break
            }
        } # End of while lopp

        precise_tail_end <- precise_tail_start + tail_end * window_size

        real_tail_start <- NA
        real_tail_end <- NA

        #plot(truncated_data, type = "l", col = "black")
        if (!is.na(precise_tail_start) & !is.na(precise_tail_end)) {
            if (read_type[seg] == 'polyA') {
                real_tail_start <- len_data - precise_tail_end + start_index
                real_tail_end <- len_data - precise_tail_start + start_index
                truncated_data <- rev(truncated_data)

                tmp <- precise_tail_start
                precise_tail_start <- len_data - precise_tail_end
                precise_tail_end <- len_data - tmp
            } else {
                real_tail_start <- start_index + precise_tail_start
                real_tail_end <- precise_tail_end + start_index
            }

            # plotting
            # plot(truncated_data, type = "l", col = "red")
            # lines(precise_tail_start:precise_tail_end,
            #       truncated_data[precise_tail_start:precise_tail_end],
            #       type = "l",
            #       col = "blue")

            # Normalize the tail length
            tail_length <- (precise_tail_end - precise_tail_start) /
                read_data$samples_per_nt

            row_num <- rca_data_tmp$id[seg]
            rca_data$tail_start[row_num] <- real_tail_start
            rca_data$tail_end[row_num] <- real_tail_end
            rca_data$tail_length[row_num] <- tail_length
        } else {
            row_num <- rca_data_tmp$id[seg]
            rca_data$tail_start[row_num] <- NA
            rca_data$tail_end[row_num] <- NA
            rca_data$tail_length[row_num] <- NA
        }
        print('s')
    } # End of for loop

    # Prepare the data
    cluster <- tail_length.x <- tail_length.y <- NULL
    rca_data_mean <- rca_data %>%
        dplyr::filter(read_type == 'polyA' | read_type == 'polyT') %>%
        dplyr::group_by(read_type, cluster) %>%
        dplyr::summarise(cluster_n = sum(!is.na(tail_length)),
                         tail_length = mean(tail_length, na.rm = TRUE)) %>%
        dplyr::mutate(file_path = file_path) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(read_type = dplyr::if_else(
            read_type == 'polyA', 'consensus_polyA', 'consensus_polyT')) %>%
        dplyr::mutate(tail_length = ifelse(is.nan(tail_length), NA, tail_length)
        )
    result <- dplyr::full_join(rca_data, rca_data_mean, by = c('read_type', 'cluster')) %>%
        dplyr::mutate(tail_length = dplyr::coalesce(tail_length.x, tail_length.y)) %>%
        dplyr::select(-tail_length.x, -tail_length.y, -id) %>%
        dplyr::mutate(samples_per_nt = samples_per_nt)
    return(result)
}
