#' Title
#'
#' @param event_data
#' @param fastq_base_number
#' @param read_type
#'
#' @return
#'
#' @examples
find_sample_index_for_fastq_base <- function(event_data, fastq_base_number, read_type){
    # get the index of the event table row containing the polyA or polyT boundart
    fqbn <- fastq_base_number - 5

    # make a cummulative sum vector
    cumm_sum <- cumsum(event_data$move)
    if (read_type=='polyT') {
        # for polyT reads take the first move that is
        # equal to the fastq_base_number
        row_index <- min(which(cumm_sum==fqbn+1))
        if (row_index==-Inf | row_index==Inf) {
            row_index <- max(which(cumm_sum==fqbn))
        }
    } else {
        # for polyA reads take the last move that is
        # equal to the fastq_base_number
        row_index <- min(which(cumm_sum==fastq_base_number-1))
        if (row_index==-Inf | row_index==Inf) {
            row_index <- min(which(cumm_sum==fastq_base_number-2))
        }
    }
    sample_index <- event_data$start[row_index]
    event_data <- cbind(event_data, cumm_sum)
    return(sample_index)
}
