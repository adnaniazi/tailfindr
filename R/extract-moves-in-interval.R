extract_moves_in_interval <- function(event_data, start, end) {
    # get the row index of the start and end raw samples in the event_data
    row_index_start <- which.min(abs(event_data$start - start))
    row_index_end <- which.min(abs(event_data$start - end))

    i <- row_index_start
    fastq_bases <- 0

    # make a fastq sequence at the most 28 bases long
    while (i <= row_index_end) {
        if (event_data$move[i] == 1) {
            fastq_bases <- fastq_bases + 1
        } else if (event_data$move[i] == 2) {
            fastq_bases <- fastq_bases + 2
        }
        i <- i + 1
    }
    return(fastq_bases)
}
