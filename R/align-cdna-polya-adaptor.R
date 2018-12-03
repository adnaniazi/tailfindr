align_cdna_polya_adaptor <- function(event_data, pri_poly_a_end, poly_a_adaptor="GAAGATAGAGCGACAGGCAAGT"){
    # get the row index of the end point of primary poly(A) tail
    row_index <- which.min(abs(event_data$start - pri_poly_a_end))

    i <- row_index
    num_bases <- 0
    fastq_bases <- ''

    # get the fastq sequence adjacent to the poly(A) tail
    # take 2 bases extra as well
    adaptor_length <- nchar(poly_a_adaptor) + 2

    while ((i < length(event_data$move)) & (num_bases < adaptor_length)){
        if (num_bases==0){
            fastq_bases <- event_data$model_state[i]
            num_bases <- 5
        } else if (event_data$move[i]==1) {
            fastq_bases <- paste(fastq_bases, substr(event_data$model_state[i], 5, 5), sep='')
            num_bases <- num_bases + 1
        } else if (event_data$move[i]==2) {
            fastq_bases <- paste(fastq_bases, substr(event_data$model_state[i], 4, 5), sep='')
            num_bases <- num_bases + 2
        }
        i <- i + 1
    }

    if (num_bases == adaptor_length+1) {
        fastq_bases <- substr(fastq_bases, 1, adaptor_length)
    }

    submat <- Biostrings::nucleotideSubstitutionMatrix(match = 1,
                                                       mismatch = -1,
                                                       baseOnly = TRUE)
    tail_adaptor_aln_score <- Biostrings::pairwiseAlignment(pattern=fastq_bases,
                                                            subject=poly_a_adaptor,
                                                            substitutionMatrix = submat,
                                                            type='local',
                                                            scoreOnly = TRUE)

    # if 6 or more bases match, then we found a good poly(A)tail
    tail_adaptor_seq <- fastq_bases
    if (tail_adaptor_aln_score > 7){
        has_valid_poly_a_tail <- TRUE
    } else {
        has_valid_poly_a_tail <- FALSE
    }
    return(list(tail_adaptor_seq = tail_adaptor_seq,
                tail_adaptor_aln_score = tail_adaptor_aln_score,
                has_valid_poly_a_tail = has_valid_poly_a_tail))
}
