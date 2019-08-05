#' Search for OligoT and reverse complement oligoT
#'
#'
#' @param fastq
#' @param ot
#' @param rc_ot
#' @param gfp
#' @param rc_gfp
#' @param ssw_before
#' @param rc_ssw_before
#' @param ssw_after
#' @param rc_ssw_after
#' @param submat
#'
#' @importFrom magrittr "%>%"
#'
#' @return a dataframe containing `start`, `stop` index of oligo in the FastQ
#' and `what_is_it` that was found (`ot` or `rc_or`). A NULL is returned if
#' no oocurance of any oligo is found.
#'
#' @export
#'
rca_search_oligo <- function(fastq,
                             ot, rc_ot,
                             gfp, rc_gfp,
                             ssw_before, rc_ssw_before,
                             ssw_after, rc_ssw_after,
                             submat) {

    # Defining thresholds
    ot_threshold <- 0.6
    rc_ot_threshold <- 0.6
    gfp_threshold <- 0.6
    rc_gfp_threshold <- 0.6
    ssw_before_threshold <-
        rc_ssw_before_threshold <-
        ssw_after_threshold <-
        rc_ssw_after_threshold <- 0.6

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

        data_list[[i]] <- list(start = start, stop = stop, what_is_it = 'OT')
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

        data_list[[i]] <- list(start = start, stop = stop, what_is_it = 'rcOT')
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


    ##########################
    ##### 3. GFP search ######
    ##########################
    fastq_tmp <- fastq

    while (TRUE) {
        searchfor <- gfp
        as <- Biostrings::pairwiseAlignment(pattern = searchfor,
                                            subject = fastq_tmp,
                                            substitutionMatrix = submat,
                                            type = type,
                                            #scoreOnly = FALSE,
                                            gapOpening = gapOpening,
                                            gapExtension = gapExtension)

        as
        nas <- as@score/searchfor@length

        if (nas < gfp_threshold) {
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

        data_list[[i]] <- list(start = start, stop = stop, what_is_it = 'GFP')
        i <- i + 1
    }

    #############################
    ##### 4. RC GFP search ######
    ############################
    fastq_tmp <- fastq
    while (TRUE) {
        searchfor <- rc_gfp
        as <- Biostrings::pairwiseAlignment(pattern = searchfor,
                                            subject = fastq_tmp,
                                            substitutionMatrix = submat,
                                            type = type,
                                            #scoreOnly = FALSE,
                                            gapOpening = gapOpening,
                                            gapExtension = gapExtension)

        as
        nas <- as@score/searchfor@length

        if (nas < rc_gfp_threshold) {
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

        data_list[[i]] <- list(start = start, stop = stop, what_is_it = 'rcGFP')
        i <- i + 1
    }


    ####################################
    ##### 5. Strand Switch Before ######
    ####################################
    fastq_tmp <- fastq
    while (TRUE) {
        searchfor <- ssw_before
        as <- Biostrings::pairwiseAlignment(pattern = searchfor,
                                            subject = fastq_tmp,
                                            substitutionMatrix = submat,
                                            type = type,
                                            #scoreOnly = FALSE,
                                            gapOpening = gapOpening,
                                            gapExtension = gapExtension)

        as
        nas <- as@score/searchfor@length

        if (nas < ssw_before_threshold) {
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

        data_list[[i]] <- list(start = start, stop = stop, what_is_it = 'SSWB')
        i <- i + 1
    }


    #######################################
    ##### 6. RC Strand Switch Before ######
    #######################################
    fastq_tmp <- fastq
    while (TRUE) {
        searchfor <- rc_ssw_before
        as <- Biostrings::pairwiseAlignment(pattern = searchfor,
                                            subject = fastq_tmp,
                                            substitutionMatrix = submat,
                                            type = type,
                                            #scoreOnly = FALSE,
                                            gapOpening = gapOpening,
                                            gapExtension = gapExtension)

        as
        nas <- as@score/searchfor@length

        if (nas < rc_ssw_before_threshold) {
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

        data_list[[i]] <- list(start = start, stop = stop, what_is_it = 'rcSSWB')
        i <- i + 1
    }


    ####################################
    ##### 7. Strand Switch After ######
    ####################################
    fastq_tmp <- fastq
    while (TRUE) {
        searchfor <- ssw_after
        as <- Biostrings::pairwiseAlignment(pattern = searchfor,
                                            subject = fastq_tmp,
                                            substitutionMatrix = submat,
                                            type = type,
                                            #scoreOnly = FALSE,
                                            gapOpening = gapOpening,
                                            gapExtension = gapExtension)

        as
        nas <- as@score/searchfor@length

        if (nas < ssw_after_threshold) {
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

        data_list[[i]] <- list(start = start, stop = stop, what_is_it = 'SSWA')
        i <- i + 1
    }


    #######################################
    ##### 8. RC Strand Switch Aftter ######
    #######################################
    fastq_tmp <- fastq
    while (TRUE) {
        searchfor <- rc_ssw_after
        as <- Biostrings::pairwiseAlignment(pattern = searchfor,
                                            subject = fastq_tmp,
                                            substitutionMatrix = submat,
                                            type = type,
                                            #scoreOnly = FALSE,
                                            gapOpening = gapOpening,
                                            gapExtension = gapExtension)

        as
        nas <- as@score/searchfor@length

        if (nas < rc_ssw_after_threshold) {
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

        data_list[[i]] <- list(start = start, stop = stop, what_is_it = 'rcSSWA')
        i <- i + 1
    }

    if (i == 1) {
        data_list <- NULL
    }

    return(data_list)
}
