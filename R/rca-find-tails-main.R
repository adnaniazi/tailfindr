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
    fastq_biostring <- Biostrings::DNAString(fastq)

    # Define the oligos adjacent to the tails
    ot <- Biostrings::DNAString("GCCTGTCGCTCTATCTTC")
    rc_ot <- Biostrings::DNAString("GAAGATAGAGCGACAGGC")

    # STEP1: Find ot and rc_ot oligos
    ot_rc_ot_df <- rca_search_oligo(fastq, ot, rc_ot, submat)

    # STEP2: Find poly(A) and poly(T) segments within the read
    polyat_df <- rca_locate_polya_polyt_segments(ot_rc_ot_df, fastq_biostring)

    # STEP3: Cluster reads
    data <- rca_cluster_read_segments(polyat_df)

    # STEP4: Create consensus of each cluster
    rca_data <- rca_create_consensus_sequence(data)

    # STEP 5: find tail lengths
    result <- rca_find_tails_per_read(rca_data = rca_data,
                                      read_data = read_data,
                                      file_path = file_path)

    return(result)
}
