#' Search for OligoT and reverse complement oligoT
#'
#' @param fastq a string. The fastQ sequence of the read.
#' @param ot a DNA string. Oligo T sequence (5' to 3')
#' @param rc_ot a DNA string. Reverese complement Oligo T sequence (5' to 3')
#' @param submat A Biostring substitution matrix
#'
#' @importFrom magrittr "%>%"
#'
#' @return a dataframe containing `start`, `stop` index of oligo in the FastQ
#' and `what_is_it` that was found (`ot` or `rc_or`). A NULL is returned if
#' no oocurance of any oligo is found.
#'
#' @export
#'
rca_search_oligo <- function(fastq, ot, rc_ot, submat) {

    # Defining thresholds
    ot_threshold <- 0.6
    rc_ot_threshold <- 0.6

    # Define parameters of alignment
    type <- 'local'
    gapOpening <- 0
    gapExtension <- 1

    ###########################
    #### 1. Oligo T search ####
    ###########################

    data_list <- list()
    fastq_tmp <- fastq

    i <- 1
    while (TRUE) {
        searchfor <- ot
        as <- Biostrings::pairwiseAlignment(pattern = searchfor,
                                            subject = fastq_tmp,
                                            substitutionMatrix = submat,
                                            type = type,
                                            gapOpening = gapOpening,
                                            gapExtension = gapExtension)

        nas <- as@score/searchfor@length

        if (nas < ot_threshold) {
            break
        }

        # substitute NNNs where alignment was found
        aligned_portion <- gsub('-', '', as@subject)
        aligned_portion_length <- nchar(aligned_portion)
        start <- as@subject@range@start
        stop <- start + aligned_portion_length
        substr(fastq_tmp,
               start = start,
               stop = stop) <-
            strrep('N', aligned_portion_length)

        data_list[[i]] <- list(start = start, stop = stop, what_is_it = 'ot')
        i <- i + 1
    }

    ##############################
    #### 2. RC Oligo T search ####
    ##############################
    fastq_tmp <- fastq
    while (TRUE) {
        searchfor <- rc_ot
        as <- Biostrings::pairwiseAlignment(pattern = searchfor,
                                            subject = fastq_tmp,
                                            substitutionMatrix = submat,
                                            type = type,
                                            gapOpening = gapOpening,
                                            gapExtension = gapExtension)

        nas <- as@score/searchfor@length

        if (nas < rc_ot_threshold) {
            break
        }

        # substitute NNNs where alignment was found
        aligned_portion <- gsub('-', '', as@subject)
        aligned_portion_length <- nchar(aligned_portion)
        start <- as@subject@range@start
        stop <- start + aligned_portion_length
        substr(fastq_tmp,
               start = start,
               stop = stop) <-
            strrep('N', aligned_portion_length)

        data_list[[i]] <- list(start = start, stop = stop, what_is_it = 'rc_ot')
        i <- i + 1
    }

    # Important to check if data_list is not empty because it will be empty
    # for erroneously fused reads
    chunk <- NULL # for CRAN
    if (length(data_list) > 0) {
        result <- data_list %>%
            purrr::map(function(.x) tibble::as_tibble(.x)) %>%
            dplyr::bind_rows(.id = "chunk") %>%
            dplyr::select(-chunk) %>%
            dplyr::arrange(start, stop)
    } else {
        result <- NULL
    }

    return(result)
}
