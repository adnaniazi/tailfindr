
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
                                         save_dir=NA,
                                         save_plots = FALSE,
                                         show_plots = FALSE,
                                         num_cores = 1,
                                         plotting_library = 'ggplot2',
                                         plot_debug = FALSE){

    if (!is.na(save_dir)) {
        dir.create(file.path(save_dir, 'plots', fsep = .Platform$file.sep))
    }

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
                                                                       show_plots = show_plots,
                                                                       save_plots = save_plots,
                                                                       save_dir = save_dir,
                                                                       plotting_library = plotting_library,
                                                                       plot_debug = plot_debug)
                                     },
                                     error=function(e){
                                         ls <- list(read_id = NA,
                                                    polya_start = NA,
                                                    polya_end = NA,
                                                    tail_length_nt = NA,
                                                    samples_per_nt = NA,
                                                    polya_fastq = NA,
                                                    file_path = file_path)
                                     })
                                 }
    close(pb)
    parallel::stopCluster(cl)

    data_list <- data.frame(data_list, stringsAsFactors = FALSE)
    row.names(data_list) <- NULL
    return(data_list)
}
