find_dna_tail <- function(fast5_dir,
                          data='cdna',
                          save_plots=FALSE,
                          show_plots=FALSE,
                          plot_debug=FALSE,
                          save_dir='~',
                          plotting_library='rbokeh'){

    message('\t- Searching for all Fast5 files...\r')
    fast5_files_list <- list.files(path=fast5_dir,
                                   pattern="\\.fast5$",
                                   recursive=TRUE,
                                   full.names=TRUE)

    message('\t- Starting a parallel cluster...\r')
    # Initiate cluster
    cl <- parallel::makeCluster(num_cores, outfile='')
    doSNOW::registerDoSNOW(cl)
    `%dopar%` <- foreach::`%dopar%`

    # progress bar
    pb <- txtProgressBar(min=1, max=length(fast5_files_list), style=3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)
    message('\t  Done!')

    #loop
    message('\t- Searching for Poly(A) and Poly(T) tails...\r')
    mcoptions <- list(preschedule=FALSE, set.seed=FALSE, cleanup=TRUE)
    data_list <-foreach::foreach(file_path = fast5_files_list,
                                 .combine = 'rbind',
                                 .options.snow = opts,
                                 .options.multicore = mcoptions) %dopar% {
                                     tryCatch({
                                         find_dna_tail_per_read(file_path=file_path,
                                                                data=data,
                                                                save_plots=save_plots,
                                                                show_plots=show_plots,
                                                                plot_debug=plot_debug,
                                                                save_dir=save_dir,
                                                                plotting_library=plotting_library)
                                     },
                                     error=function(e){
                                         ls <- list(read_id = NA,
                                                    read_type = NA,
                                                    tail_is_valid = NA,
                                                    tail_start = NA,
                                                    tail_end = NA,
                                                    samples_per_nt = NA,
                                                    tail_length = NA,
                                                    file_path = file_path,
                                                    has_precise_boundary = NA)
                                     })
                                 }
    close(pb)
    parallel::stopCluster(cl)

    data_list <- data.frame(data_list, stringsAsFactors = FALSE)
    data_list$read_id <- as.character(data_list$read_id)
    return(data_list)
}
