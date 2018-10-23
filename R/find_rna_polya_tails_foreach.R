
#' Title
#'
#' @param fast5_files_list
#' @param save_dir
#' @param csv_file_name
#' @param save_plots
#' @param show_plots
#' @param num_cores
#' @param plotting_library
#' @param plot_debug
#'
#' @return
#' @export
#'
#' @examples
find_rna_polya_tails_foreach <- function(fast5_files_list,
                                         save_dir,
                                         save_plots=FALSE,
                                         show_plots=FALSE,
                                         num_cores=1,
                                         plotting_library='ggplot2',
                                         plot_debug=FALSE){

    message('\t- Starting a parallel cluster...\r')
    # Initiate cluster
    cl <- parallel::makeCluster(num_cores)
    doSNOW::registerDoSNOW(cl)
    `%dopar%` <- foreach::`%dopar%`
    `%do%` <- foreach::`%do%`

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
                                         find_rna_polya_tail_per_read(file_path,
                                                                       show_plots=show_plots,
                                                                       save_plots=save_plots,
                                                                       save_dir=save_dir,
                                                                       plotting_library=plotting_library,
                                                                       plot_debug=plot_debug)
                                     },
                                     error=function(e){
                                         ls <- list(read_id = NA,
                                                    poly_a_start = NA,
                                                    poly_a_end = NA,
                                                    poly_a_fastq = NA,
                                                    poly_a_length_in_nucleotides_1 = NA,
                                                    poly_a_length_in_nucleotides_2 = NA,

                                                    non_poly_a_seq_start = NA,
                                                    non_poly_a_seq_end = NA,
                                                    moves_in_non_poly_a_region = NA,

                                                    sampling_rate = NA,
                                                    file_path=file_path)
                                     })
                                 }
    close(pb)
    parallel::stopCluster(cl)

    data_list <- data.frame(data_list)
    data_list$read_id <- as.character(data_list$read_id)
    return(data_list)
}
