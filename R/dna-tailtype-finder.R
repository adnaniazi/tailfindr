#' Title
#'
#' @param file_path
#'
#' @return
#' @export
#'
#' @examples
dna_tailtype_finder <- function(file_path) {

    read_data <- extract_read_data_hdf5r(file_path, plot_debug=FALSE)
    fastq <- Biostrings::DNAString(read_data$fastq)

    fa <- Biostrings::DNAString('GGCGTCTGCTTGGGTGTTTAACCTTTTTTTTTTAATGTACTTCGTTCAGTTACGTATTGCT')
    ea <- Biostrings::DNAString('GCAATACGTAACTGAACGAAGT')

    rc_fa <- Biostrings::reverseComplement(fa)
    rc_ea <- Biostrings::reverseComplement(ea)

    match = 1
    mismatch = -1
    type='local'
    gapOpening=0
    gapExtension=1

    submat <- Biostrings::nucleotideSubstitutionMatrix(match = match,
                                                       mismatch = mismatch,
                                                       baseOnly = TRUE)

    aln_score_gfp <- Biostrings::pairwiseAlignment(pattern=fa,
                                                   subject=fastq,
                                                   substitutionMatrix = submat,
                                                   type=type,
                                                   scoreOnly = FALSE,
                                                   gapOpening=gapOpening,
                                                   gapExtension=gapExtension)


    # get the start chunk
    aln_score_gfp
    aln_score_gfp

}

