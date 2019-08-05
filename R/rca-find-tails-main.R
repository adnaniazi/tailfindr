#' Title
#'
#' @param file_path
#' @param read_id_fast5_file
#' @param basecall_group
#' @param save_plots
#' @param show_plots
#' @param plot_debug
#' @param save_dir
#' @param plotting_library
#' @param multifast5
#' @param basecalled_with
#' @param model
#' @param ...
#'
#' @return
#' @export
#'
rca_find_tails_main <- function(file_path,
                                read_id_fast5_file = NA,
                                basecall_group = 'Basecall_1D_000',
                                save_plots = save_plots,
                                show_plots = show_plots,
                                plot_debug = FALSE,
                                save_dir = save_dir,
                                plotting_library = plotting_library,
                                multifast5 = multifast5,
                                basecalled_with = basecalled_with,
                                model = model,
                                ...) {

    # get the substitution matrix -- if passed from outside
    if (length(list(...)) > 0) {
        lst <- list(...)
        if ("dna_opts" %in% names(lst)) {
            dna_opts <- lst$dna_opts
            match <- dna_opts$match
            mismatch <- dna_opts$mismatch
            type <- dna_opts$type
            gapOpening <- dna_opts$gapOpening
            gapExtension <- dna_opts$gapExtension
            submat <- dna_opts$submat
        }
    } else {
        # otherwise make one
        match <- 1
        mismatch <- -1
        type <- 'local'
        gapOpening <- 0
        gapExtension <- 1
        submat <- Biostrings::nucleotideSubstitutionMatrix(match = match,
                                                           mismatch = mismatch)
    }
    # extract read data
    read_data <- extract_read_data(file_path,
                                   read_id_fast5_file,
                                   plot_debug,
                                   basecalled_with,
                                   basecall_group,
                                   multifast5,
                                   model,
                                   plotting_library)
    fastq <- read_data$fastq
    fastq_length <- nchar(fastq)
    #fastq_biostring <- Biostrings::DNAString(fastq)

    # Define the oligos adjacent to the tails
    ot <- Biostrings::DNAString("GCCTGTCGCTCTATCTTC")
    rc_ot <- Biostrings::DNAString("GAAGATAGAGCGACAGGC")
    gfp <- Biostrings::DNAString('GGATCACTCTCGGCATGGACGAGCTGTACAAGTAG')
    rc_gfp <- Biostrings::reverseComplement(gfp)
    ssw_before <- Biostrings::DNAString("CTGTTGGTGCTGATATTGC")
    rc_ssw_before <-  Biostrings::reverseComplement(ssw_before)
    ssw_after <- Biostrings::DNAString("CGATGTACGGCCGGG")
    rc_ssw_after <-  Biostrings::reverseComplement(ssw_after)

    # STEP1: Find ot and rc_ot oligos
    data_list <- rca_search_oligo(fastq,
                                  ot, rc_ot,
                                  gfp, rc_gfp,
                                  ssw_before, rc_ssw_before,
                                  ssw_after, rc_ssw_after, submat)


    # TODO: IMPORTANT do this only if data_list is not empty
    chunk <- NULL
    if (!is.null(data_list)) {
        result <- data_list
        result <- purrr::map(result, function(.x) tibble::as_tibble(.x))
        result <- dplyr::bind_rows(result, .id = "chunk")
        result <- dplyr::select(result, -chunk)
        result <- result %>% dplyr::arrange(start, stop)
        content_string <- paste(result$what_is_it, collapse = '-')
        fastq_start_string <- paste(result$start, collapse = '-')

        return(list(read_id = read_data$read_id,
                    content_string = content_string,
                    fastq_start_string = fastq_start_string,
                    fastq_length = fastq_length,
                    file_path = file_path))
    } else {
        return(list(read_id = read_data$read_id,
                    content_string = NA,
                    fastq_start_string = NA,
                    fastq_length = fastq_length,
                    file_path = file_path))
    }

}
