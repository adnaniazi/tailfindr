#' Find Poly(A) tails in FAST5 reads
#'
#' This function gets a directory of forward-strand FAST5 reads, searches
#' for poly(A) tails in these reads, and then returns a dataframe that contains
#' the data about poly-(A) tails found
#'
#' @param fast5_files_list a character list.
#' @param save_dir a character string. Path of the directory in which to save the
#' results csv file and plots (optional)
#' @param csv_file_name a character string. Name of the csv file in which to store the results.
#' @param save_plots A logical [FALSE]. Set it to TRUE if you want to
#' save a plot of poly(A) tail in each read, otherwise set it to FALSE.
#' @param show_plots A logical [FALSE]. Set it to TRUE if you want to display a
#' plot of poly(A) tails in each read, otherwise set it to FALSE. Setting this
#' option to TRUE could cause the algorithm to run very slowly.
#' @param num_cores
#'
#' @return A dataframe containing information about poly(A) tail lengths of each read
#' @export
#'
#' @examples
#' df <- find_cdna_polya_tails_batch_parallel('/FORWARD/STRAND/FAST5/FILES/DIRECTORY', 'SAVE/DIR', 'polya-tail-data.csv')
find_cdna_polya_tails_foreach <- function(fast5_files_list,
                                                 save_dir,
                                                 csv_file_name,
                                                 save_plots=FALSE,
                                                 show_plots=FALSE,
                                                 num_cores=1){

    message('\t- Starting a parallel cluster...\r')
    # Initiate cluster
    cl <- parallel::makeCluster(num_cores)
    doSNOW::registerDoSNOW(cl)
    `%dopar%` <- foreach::`%dopar%`

    # progress bar
    pb <- txtProgressBar(min=1, max=length(fast5_files_list), style=3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)
    message('\t  Done!')

    #loop
    message('\t- Searching for Poly(A) tails...\r')
    mcoptions <- list(preschedule=FALSE, set.seed=FALSE, cleanup=TRUE)
    data_list <-foreach::foreach(file_path = fast5_files_list,
                                 .combine = 'rbind',
                                 .options.snow = opts,
                                 .options.multicore = mcoptions) %dopar% {
        tryCatch({
            find_cdna_polya_tail_per_read(file_path,
                                          show_plots=show_plots,
                                          save_plots=save_plots,
                                          save_dir=save_dir)
        },
        error=function(e){
            ls <- list(read_id=NA,
                       pri_poly_a_start=NA,
                       pri_poly_a_end=NA,
                       pri_poly_a_fastq=NA,
                       gap1_start=NA,
                       gap1_end=NA,
                       gap1_fastq=NA,
                       sec1_poly_a_start=NA,
                       sec1_poly_a_end=NA,
                       sec1_poly_a_fastq=NA,
                       gap2_start=NA,
                       gap2_end=NA,
                       gap2_fastq=NA,
                       sec2_poly_a_start=NA,
                       sec2_poly_a_end=NA,
                       sec2_poly_a_fastq=NA,
                       non_poly_a_seq_start = NA,
                       non_poly_a_seq_end = NA,
                       moves_in_non_poly_a_region = NA,
                       sampling_rate=NA,
                       cdna_poly_a_read_type='Fatal Error',
                       tail_adaptor=NA,
                       has_valid_poly_a_tail=FALSE,
                       file_path=file_path)
        })
    }
    close(pb)
    parallel::stopCluster(cl)

    data_list <- data.frame(data_list)
    data_list$read_id <- as.character(data_list$read_id)
    return(data_list)
}
