#' Cluster segments found in a long RCA reads
#'
#' Sometimes two or more completely different transcripts may be found in a
#' read. Multiple copies of those then need to be clustered together.
#'
#' @param polyat_df a dataframe as produced by the
#' `rca_locate_polya_polyt_segments` function
#'
#' @return
#' @export
#'
rca_cluster_read_segments <- function(polyat_df) {

    if (is.null(polyat_df)) {
        return(NULL)
    }

    read_type <- NULL # for CRAN

    # make polyA and polyT-containing dataframes
    polya_df <- dplyr::filter(polyat_df, read_type == 'polyA')
    polyt_df <- dplyr::filter(polyat_df, read_type == 'polyT')

    # check if there are polyA reads in the first place
    if (dim(polya_df)[1] != 0) {
        # make DNA string sets of polyA reads
        polya_seqs <- Biostrings::DNAStringSet(polya_df$fastq_segment)

        # Cluster polyA reads based on similarity
        polya_cluster <- DECIPHER::IdClusters(myDistMatrix = NULL,
                                              method = "inexact",
                                              myXStringSet = polya_seqs,
                                              cutoff = 0.9,
                                              processors = NULL,
                                              verbose = FALSE)

        # add the cluster column to the df
        polya_df <- dplyr::mutate(polya_df, cluster = polya_cluster$cluster)

        # Make a list of polyA read clusters,
        # each list contains all the sequences of that cluster
        polya_cluster_list <- list()
        for (i in sort(unique(polya_cluster$cluster))) {
            polya_cluster_list[[i]] <-
                polya_df$fastq_segment[polya_cluster$cluster == i]
        }
    } else {
        polya_cluster_list <- NULL
    }

    # check if there are polyT reads in the first place
    if (dim(polyt_df)[1] != 0) {
        # make DNA string sets of polyA reads
        polyt_seqs <- Biostrings::DNAStringSet(polyt_df$fastq_segment)

        polyt_cluster <- DECIPHER::IdClusters(myDistMatrix = NULL,
                                              method = "inexact",
                                              myXStringSet = polyt_seqs,
                                              cutoff = 0.9,
                                              processors = NULL,
                                              verbose = FALSE)

        # add the cluster column to the df
        polyt_df <- dplyr::mutate(polyt_df, cluster = polyt_cluster$cluster)

        # Make a list of polyA read clusters,
        # each list contains all the sequences of that cluster
        polyt_cluster_list <- list()
        for (i in sort(unique(polyt_cluster$cluster))) {
            polyt_cluster_list[[i]] <-
                polyt_df$fastq_segment[polyt_cluster$cluster == i]
        }
    } else {
        polyt_cluster_list <- NULL
    }

    # combine the polyA and polyT data frames into one
    polyat_df <- rbind(polyt_df, polya_df)

    return(
        list(
            polya_cluster_list = polya_cluster_list,
            polyt_cluster_list = polyt_cluster_list,
            polyat_df = polyat_df
            )
    )
}
