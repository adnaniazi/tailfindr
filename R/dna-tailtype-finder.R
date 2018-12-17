#' Title
#'
#' @param file_path
#' @param plot_debug
#'
#' @return
#' @export
#'
#' @examples
dna_tailtype_finder <- function(file_path, data='cdna', plot_debug=FALSE) {

    match = 1
    mismatch = -1
    type='local'
    gapOpening=0
    gapExtension=1
    submat <- Biostrings::nucleotideSubstitutionMatrix(match = match,
                                                       mismatch = mismatch,
                                                       baseOnly = TRUE)

    read_data <- extract_read_data_hdf5r(file_path, plot_debug=plot_debug)
    event_data <- read_data$event_data
    fastq <- read_data$fastq

    fa <- Biostrings::DNAString('GGCGTCTGCTTGGGTGTTTAACCTTTTTTTTTTAATGTACTTCGTTCAGTTACGTATTGCT')
    ea <- Biostrings::DNAString('GCAATACGTAACTGAACGAAGT')
    # CDNA
    if (data == 'cdna') {
        fp <- Biostrings::DNAString('TTTCTGTTGGTGCTGATATTGCTGCCATTACGGCCGGG')
        ep <- Biostrings::DNAString('GAAGATAGAGCGACAGGCAAGT')
    } else if  (data == 'pcr-dna') {
        fp <- Biostrings::DNAString('ATTTAGGTGACACTATAGCGCTCCATGCAAACCTGTC')
        ep <- Biostrings::DNAString('CGTTGCCGCCCGGACTC')
    }

    rc_fp <- Biostrings::reverseComplement(fp)
    rc_ep <- Biostrings::reverseComplement(ep)

    as_fp <- Biostrings::pairwiseAlignment(pattern=fp,
                                           subject=Biostrings::DNAString(substr(fastq, start=1, stop=min(100, nchar(fastq)))),
                                           substitutionMatrix = submat,
                                           type=type,
                                           scoreOnly = FALSE,
                                           gapOpening=gapOpening,
                                           gapExtension=gapExtension)

    as_rc_ep <- Biostrings::pairwiseAlignment(pattern=rc_ep,
                                              subject=Biostrings::DNAString(substr(fastq, start=1, stop=min(70, nchar(fastq)))),
                                              substitutionMatrix = submat,
                                              type=type,
                                              scoreOnly=FALSE,
                                              gapOpening=gapOpening,
                                              gapExtension=gapExtension)


    nas_fp <- as_fp@score/fp@length
    nas_rc_ep <- as_rc_ep@score/rc_ep@length
    nas_ep <- NA # for max remove later
    nas_rc_fp <- NA # for max reomve later
    has_precise_boundary <- FALSE
    bases_to_match <- 3

    # check the front primer and rev comp end primer score to decide
    # between polyA and polyT reads.
    if (nas_fp > nas_rc_ep & nas_fp > 0.6){
        read_type <- 'polyA'
        # check if there is the end primer at the end of the polyA read
        # adjacent to the polyA tail
        as_ep <- Biostrings::pairwiseAlignment(pattern=ep,
                                               subject=Biostrings::DNAString(substr(fastq, start=max(nchar(fastq)-50, 0), stop=nchar(fastq))),
                                               substitutionMatrix = submat,
                                               type=type,
                                               scoreOnly=FALSE,
                                               gapOpening=gapOpening,
                                               gapExtension=gapExtension)

        nas_ep <- as_ep@score/ep@length
        tail_is_valid <- ifelse(nas_ep>0.6, T, F)
        if (tail_is_valid) {
            polya_end_fastq <- as_ep@subject@range@start + nchar(fastq)-50
            polya_rough_end <- find_sample_index_for_fastq_base(read_data$event_data, polya_end_fastq, read_type)
            # for max remove later
            as_rc_fp <- Biostrings::pairwiseAlignment(pattern=rc_fp,
                                                   subject=Biostrings::DNAString(substr(fastq, start=max(nchar(fastq)-50, 0), stop=nchar(fastq))),
                                                   substitutionMatrix = submat,
                                                   type=type,
                                                   scoreOnly=FALSE,
                                                   gapOpening=gapOpening,
                                                   gapExtension=gapExtension)

            # check if we have captured the end of the tail perfectly
            # by finding out if we captured the bases adjacent to the tails
            if (substr(as_ep@subject,
                       start=1,
                       stop=bases_to_match) == substr(ep,
                                                        start=1,
                                                        stop=bases_to_match)) {
                has_precise_boundary <- TRUE
            }

            # for max remove later
            nas_rc_fp <- as_rc_fp@score/rc_fp@length

        } else {
            polya_end_fastq <- NA
            polya_rough_end <- NA
        }
        polyt_start_fastq <- NA
        polyt_rough_start <- NA

    # Check if it is a PolyT tail
    } else if (nas_fp < nas_rc_ep & nas_rc_ep > 0.6) {
        read_type <- 'polyT'
        tail_is_valid <- T
        polyt_start_fastq <- as_rc_ep@subject@range@start + as_rc_ep@subject@range@width
        polyt_rough_start <- find_sample_index_for_fastq_base(read_data$event_data, polyt_start_fastq, read_type)
        polya_end_fastq <- NA
        polya_rough_end <- NA

        # check if we have captured the end of the tail perfectly
        # by finding out if we captured the bases adjacent to the tails
        if (substr(as_rc_ep@subject,
                   start=nchar(as.character(as_rc_ep@subject))-bases_to_match+1,
                   stop=nchar(as.character(as_rc_ep@subject))) == substr(rc_ep,
                                                                         start=rc_ep@length-bases_to_match+1,
                                                                         stop=rc_ep@length)) {
            has_precise_boundary <- TRUE
        }


    # if the above two checks fail then it is an invalid read
    } else {
        read_type <- 'invalid'
        tail_is_valid <- F
        polyt_start_fastq <- NA
        polyt_rough_start <- NA
        polya_end_fastq <- NA
        polya_rough_end <- NA
        has_precise_boundary <- NA
    }

    df <- data.frame(read_data$raw_data)
    # p <- ggplot2::ggplot(data=df) +
    #     ggplot2::geom_line(ggplot2::aes(x=c(1:length(read_data$raw_data)), y=read_data.raw_data))
    #
    # if (!is.na(polyt_rough_start)){
    #     p <- p + ggplot2::geom_vline(xintercept = polyt_rough_start, color='red')
    # }
    # if (!is.na(polya_rough_end)){
    #     p <- p + ggplot2::geom_vline(xintercept = polya_rough_end, color='green')
    # }
    #
    # print(p)

    # nchar(fastq)
    return(list(read_data=read_data, # for max comment it out
                read_type=read_type,
                tail_is_valid=tail_is_valid,
                polya_end_fastq=polya_end_fastq,
                polyt_start_fastq=polyt_start_fastq,
                polya_end=polya_rough_end,
                polyt_start=polyt_rough_start,
                has_precise_boundary=has_precise_boundary,
                nas_fp=nas_fp,
                nas_rc_ep=nas_rc_ep,
                nas_ep=nas_ep,
                nas_rc_fp=nas_rc_fp))
}

