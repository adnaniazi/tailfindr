#' Title
#'
#' @param event_data
#' @param fastq_start_index
#' @param fastq_end_index
#' @param read_type
#' @param raw_data_len
#' @param samples_per_nt
#'
#' @return
#' @export
#'
#' @examples
rca_get_segment_boundaries_in_samples <- function(event_data,
                                                 fastq_start_index,
                                                 fastq_end_index,
                                                 read_type,
                                                 raw_data_len,
                                                 samples_per_nt) {

    # some padding before the start of one of the tail ends near ot
    # and rc_ot
    extra_bases_near_oligos <- 1
    search_window_size_nt <- 400

    if (read_type == 'polyA') {
        # get end sample index
        sample_end_idx <-
            find_sample_index_for_fastq_base(
                event_data = event_data,
                fastq_base_number = fastq_end_index,
                read_type = 'polyA'
            )

        # shift end sample index X nucleotides for safety
        sample_end_idx <- floor(
            min(raw_data_len,
                sample_end_idx + extra_bases_near_oligos * samples_per_nt)
        )

        # get start sample index (go to the left of poly(A) tail by 600 nt)
        sample_start_idx <- floor(
            max(1, sample_end_idx - search_window_size_nt * samples_per_nt)
        )

    } else if (read_type == 'polyT') {
        sample_start_idx <-
            find_sample_index_for_fastq_base(
                event_data = event_data,
                fastq_base_number = fastq_start_index,
                read_type = 'polyT'
            )

        # shift sample_start_idx by X nucleotides for safety
        sample_start_idx <- floor(
            max(1,
                sample_start_idx - extra_bases_near_oligos * samples_per_nt)
        )

        # get end sample index (go fto the right of poly(T) tail by 600 nt)
        sample_end_idx <- floor(
            min(raw_data_len,
                sample_start_idx + search_window_size_nt * samples_per_nt)
        )
    }

    return(list(
        sample_start_idx = sample_start_idx,
        sample_end_idx = sample_end_idx
    ))

}
