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
    aln_score <- Biostrings::pairwiseAlignment(pattern=fastq_bases,
                                               subject=poly_a_adaptor,
                                               substitutionMatrix = submat,
                                               type='local',
                                               scoreOnly = TRUE)

    # if 6 or more bases match, then we found a good poly(A)tail
    if (aln_score > 5){
        tail_adaptor <- paste('Tail poly_a_adaptor found; aln score: ', aln_score, '; poly_a_adaptor seq: ', fastq_bases, sep='')
        has_valid_poly_a_tail <- TRUE
    } else {
        tail_adaptor <- paste('Tail poly_a_adaptor absent; aln score: ', aln_score, '; poly_a_adaptor seq: ', fastq_bases, sep='')
        has_valid_poly_a_tail <- FALSE
    }
    return(c(tail_adaptor, has_valid_poly_a_tail))
}
