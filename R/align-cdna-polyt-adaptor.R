align_cdna_polyt_adaptor <- function(event_data, pri_poly_t_start, polyt_adaptor="ACTTGCCTGTCGCTCTATCTTC"){
    # get the row index of the end point of primary poly(A) tail
    row_index <- which.min(abs(event_data$start - pri_poly_t_start))

    i <- row_index
    num_bases <- 0
    fastq_bases <- ''

    # get the fastq sequence adjacent to the poly(A) tail
    # take 2 bases extra as well
    adaptor_length <- nchar(polya_adaptor) + 2

    while ((i > 0) & (num_bases < adaptor_length)){
        if (event_data$move[i]==1) {
            fastq_bases <- paste(substr(event_data$model_state[i], 5, 5), fastq_bases, sep='')
            num_bases <- num_bases + 1
        } else if (event_data$move[i]==2) {
            fastq_bases <- paste(substr(event_data$model_state[i], 4, 5), fastq_bases, sep='')
            num_bases <- num_bases + 2
        }
        i <- i - 1
    }

    if (num_bases == adaptor_length+1) {
        fastq_bases <- substr(fastq_bases, 1, adaptor_length)
    }

    submat <- Biostrings::nucleotideSubstitutionMatrix(match = 1,
                                                       mismatch = -1,
                                                       baseOnly = TRUE)
    aln_score <- Biostrings::pairwiseAlignment(pattern=fastq_bases,
                                               subject=polyt_adaptor,
                                               substitutionMatrix = submat,
                                               type='local',
                                               scoreOnly = TRUE)

    # if 6 or more bases match, then we found a good poly(A)tail
    if (aln_score > 5){
        tail_adaptor <- paste('Tail polyt_adaptor found; aln score: ', aln_score, '; polyt_adaptor seq: ', fastq_bases, sep='')
        has_valid_poly_t_tail <- TRUE
    } else {
        tail_adaptor <- paste('Tail polyt_adaptor absent; aln score: ', aln_score, '; polyt_adaptor seq: ', fastq_bases, sep='')
        has_valid_poly_t_tail <- FALSE
    }
    return(c(tail_adaptor, has_valid_poly_t_tail))
}
