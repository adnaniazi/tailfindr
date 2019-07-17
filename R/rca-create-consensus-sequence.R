#' Create consensus sequence from multiple sequence
#'
#' As RCA read contains multiple copies of a transcript isoform, this function
#' create a consensus sequence
#'
#' @param data
#'
#' @return
#' @export
#'
rca_create_consensus_sequence <- function(data){

    polyat_df <- data$polyat_df
    polya_cluster_list <- data$polya_cluster_list
    polyt_cluster_list <- data$polyt_cluster_list

    gapOpening <- 0
    gapExtension <- 1

    tool <- 'decipher'
    # msa tool has this wierd error
    # internalRsequence.aln and internalRsequence.dnd not found

    # DO the alignment of all clusters
    for (i in seq_along(polya_cluster_list)) {
        if (length(polya_cluster_list[[i]]) > 1) {
            if (tool == 'msa') {
                alignments <- msa::msa(Biostrings::DNAStringSet(polya_cluster_list[[i]]),
                                       method = 'ClustalW',
                                       gapOpening = gapOpening,
                                       gapExtension = gapExtension)
                consensus <- msa::msaConsensusSequence(alignments, ignoreGaps = FALSE)
            } else {
                alignments <- DECIPHER::AlignSeqs(Biostrings::DNAStringSet(polya_cluster_list[[i]]),
                                                  gapOpening = gapOpening,
                                                  gapExtension = gapExtension,
                                                  verbose = FALSE)
                consensus <- DECIPHER::ConsensusSequence(alignments,
                                                         ambiguity = FALSE,
                                                         threshold=0.3,
                                                         minInformation=0.6,
                                                         noConsensusChar="-")
                consensus <- as.character(consensus)
            }
            consensus <- gsub('[[:punct:]]', '', consensus)
        } else {# if only one element in cluster then clustering fails
            consensus <- as.character(
                Biostrings::DNAString(polya_cluster_list[[i]])
            )
        }
        polyat_df[nrow(polyat_df) + 1, ] = list(start = NA,
                                                end = NA,
                                                fastq_segment = consensus,
                                                read_type = 'consensus_polyA',
                                                cluster = i)
    }

    for (i in seq_along(polyt_cluster_list)) {
        if (length(polyt_cluster_list[[i]]) > 1) {
            if (tool == 'msa') {
                alignments <- msa::msa(Biostrings::DNAStringSet(polyt_cluster_list[[i]]),
                                       method = 'ClustalW',
                                       gapOpening = gapOpening,
                                       gapExtension = gapExtension)
                consensus <- msa::msaConsensusSequence(alignments, ignoreGaps = FALSE)
            } else {
                alignments <- DECIPHER::AlignSeqs(Biostrings::DNAStringSet(polyt_cluster_list[[i]]),
                                                  gapOpening = gapOpening,
                                                  gapExtension = gapExtension,
                                                  verbose = FALSE)
                consensus <- DECIPHER::ConsensusSequence(alignments,
                                                         ambiguity = FALSE,
                                                         threshold=0.3,
                                                         minInformation=0.6,
                                                         noConsensusChar="-")
                consensus <- as.character(consensus)
            }
            consensus <- gsub('[[:punct:]]', '', consensus)
        } else {# if only one element in cluster then clustering fails
            consensus <- as.character(
                Biostrings::DNAString(polyt_cluster_list[[i]])
            )
        }

        polyat_df[nrow(polyat_df) + 1, ] = list(start = NA,
                                                end = NA,
                                                fastq_segment = consensus,
                                                read_type = 'consensus_polyT',
                                                cluster = i)
    }

    return(polyat_df)
}
