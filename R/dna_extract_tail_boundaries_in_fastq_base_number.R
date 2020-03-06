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
dna_extract_tail_boundaries_in_fastq_base_number <- function(event_data, start, end, fastq) {
    # get the row index of the start and end raw samples in the event_data
    row_index_start <- which.min(abs(event_data$start - start))
    row_index_end <- which.min(abs(event_data$start - end))

    move_cummsum <- cumsum(event_data$move)
    event_data <- cbind(event_data, move_cummsum)
    start_base_index <- event_data$move_cummsum[row_index_start]
    end_base_index <- event_data$move_cummsum[row_index_end]
    #fastq_length <- nchar(fastq)
    #end_base_index <- fastq_length - end_base_index + 1
    #start_base_index <- fastq_length - start_base_index + 1

    return(list(
        start_base_index = start_base_index,
        end_base_index = end_base_index
    ))
}
