#' Title
#'
#' @param event_data
#' @param start
#' @param end
#' @param fastq
#'
#' @return
#'
#' @examples
extract_fasta_between_start_and_end_sample_indices <- function(event_data, start, end, fastq) {
    # get the row index of the start and end raw samples in the event_data
    row_index_start <- which.min(abs(event_data$start - start))
    row_index_end <- which.min(abs(event_data$start - end))

    move_cummsum <- cumsum(event_data$move)
    event_data <- cbind(event_data, move_cummsum)
    start_base_index <- event_data$move_cummsum[row_index_start]
    end_base_index <- event_data$move_cummsum[row_index_end]
    fastq_length <- nchar(fastq)
    fasta_bases <- substr(fastq, start=start_base_index, stop=end_base_index)
    return(list(fasta_bases = fasta_bases,
                start_base_index = start_base_index)
           )
}
