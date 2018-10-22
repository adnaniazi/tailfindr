extract_fastq_in_interval <- function(event_data, start, end){
    # get the row index of the start and end raw samples in the event_data
    row_index_start <- which.min(abs(event_data$start - start))
    row_index_end <- which.min(abs(event_data$start - end))

    i <- row_index_start
    fastq_bases <- ''

    # make a fastq sequence
    while (i <= row_index_end){
        if (i == row_index_start){
            fastq_bases <- event_data$model_state[i]
        } else if (event_data$move[i] == 1) {
            fastq_bases <- paste(fastq_bases, substr(event_data$model_state[i], 5, 5), sep='')
        } else if (event_data$move[i] == 2) {
            fastq_bases <- paste(fastq_bases, substr(event_data$model_state[i], 4, 5), sep='')
        }
        i <- i + 1
    }
    return(fastq_bases)
}
