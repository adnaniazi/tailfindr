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
    #moves_cumsum <- cumsum(event_data$move)
    #event_data <- cbind(event_data, moves_cumsum)
    model_state_length <- nchar(event_data$model_state[1])

    if (model_state_length > 1) {
        fqbn <- fastq_base_number - model_state_length + 1
    } else {
        fqbn <- fastq_base_number
    }

    # make a cummulative sum vector
    cumm_sum <- cumsum(event_data$move)
    if (read_type=='polyT') {
        # for polyT reads take the first move that is
        # equal to the fastq_base_number
        row_index <- min(which(cumm_sum==fqbn))
        if (row_index==-Inf | row_index==Inf) {
            row_index <- max(which(cumm_sum==fqbn-1))
        }
    } else {
        # for polyA reads take the last move that is
        # equal to the fastq_base_number
        row_index <- min(which(cumm_sum==fqbn)) - 1
        if (row_index==-Inf | row_index==Inf) {
            row_index <- max(which(cumm_sum==fqbn-1))
        }
    }
    sample_index <- event_data$start[row_index]
    return(sample_index)
}
