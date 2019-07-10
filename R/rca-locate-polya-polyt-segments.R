#' Find polyA/PolyT segments within the long RCA reads
#'
#' @param ot_rc_ot_df a dataframe as returned by the rca_search_oligo function
#' @param fastq a DNA string. The fastQ sequence of the read
#' @importFrom magrittr "%>%"
#' @return
#' @export
#'
rca_locate_polya_polyt_segments <- function(ot_rc_ot_df,
                                            fastq) {

    if (is.null(ot_rc_ot_df)) {
        return(NULL)
    }

    what_is_it <- ot_rc_ot_df$what_is_it
    starts <- ot_rc_ot_df$start
    stops <- ot_rc_ot_df$stop
    minimum_segment_length <- 50 # segment found should be larger than this

    i <- 1
    polyat_seqs <- list()
    polyat_seq_count <- 0

    # Recipe:
    # PolyT segments are those which are within two consecutive ot occurances
    # PolyA segments are those which are within two consecutive rc_ot occurances

    len_what_is_it <- length(what_is_it)
    while (i < len_what_is_it) {
        if (what_is_it[i] == 'ot' & what_is_it[i + 1] == 'ot' ) {
            # extract the RC transcript sequence in between the oligo T segments
            start_site <- stops[i]
            stop_site <- starts[i + 1]
            fastq_segment <- fastq[start_site:stop_site]

            # if it is too short a segment, then it is probably some erroneous
            # fusion event, ignor it then
            if (stop_site - start_site > minimum_segment_length) {
                polyat_seq_count <- polyat_seq_count + 1
                polyat_seqs[[polyat_seq_count]] <-
                    list(
                        start = start_site,
                        end = stop_site,
                        fastq_segment = as.character(fastq_segment),
                        read_type = 'polyT'
                    )
            }
        } else if (what_is_it[i] == 'rc_ot' & what_is_it[i + 1] == 'rc_ot' ) {
            # extract the transcript sequence in between the RC oligo T segments
            start_site <- stops[i]
            stop_site <- starts[i + 1]
            fastq_segment <- fastq[start_site:stop_site]

            # if it is too short a segment, then it is probably some erroneous
            # fusion event, ignor it then
            if (stop_site - start_site > minimum_segment_length) {
                polyat_seq_count <- polyat_seq_count + 1
                polyat_seqs[[polyat_seq_count]] <-
                    list(
                        start = start_site,
                        end = stop_site,
                        fastq_segment = as.character(fastq_segment),
                        read_type = 'polyA'
                    )
            }
        }
        i <- i + 1
    }

    # Important to check if data_list is not empty because it will be empty
    # for erroneously fused reads
    if (length(polyat_seqs) > 0) {
        result <- polyat_seqs %>%
            purrr::map(function(.x) tibble::as_tibble(.x)) %>%
            dplyr::bind_rows(.id = "chunk") %>%
            dplyr::select(-chunk) %>%
            dplyr::arrange(start, end)
    } else {
        result <- NULL
    }

    return(result)
}

