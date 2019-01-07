#' Title
#'
#' @param fast5_dir
#' @param data
#' @param save_dir
#' @param csv_filename
#' @param num_cores
#' @param save_plots
#' @param show_plots
#' @param plot_debug
#' @param plotting_library
#'
#' @return
#' @export
#'
#' @examples
find_dna_tails <- function(fast5_dir,
                           data='cdna',
                           save_dir='~',
                           csv_filename='dna_tail.csv',
                           num_cores=1,
                           save_plots=FALSE,
                           show_plots=FALSE,
                           plot_debug=FALSE,
                           plotting_library='rbokeh'){

    # Try to create the save directory
    if (!dir.exists(file.path(save_dir))) {
        message('\t- Save dir does not exist. Trying to create it...\r')
        tryCatch({
            dir.create(file.path(save_dir, fsep = .Platform$file.sep))
            message('\t  Done!')
        },
        error=function(e){
            message('\t  Failed to create the save dir. Results will be stored in "~/" directory instead\r')
            save_dir = '~/'
        })
    }

    # Create a sub-direcotry to save all the plots
    if (save_plots){
        message('\t- Creating a sub-directory to save the plots in...\r')
        dir.create(file.path(save_dir, 'plots', fsep = .Platform$file.sep))
        message('\t  Done! All plots will be saved in the following direcotry:\r')
        message(file.path(save_dir, 'plots', fsep = .Platform$file.sep))
    }

    message('\t- Searching for all Fast5 files...\r')
    fast5_files_list <- list.files(path=fast5_dir,
                                   pattern="\\.fast5$",
                                   recursive=TRUE,
                                   full.names=TRUE)
    num_files <- length(fast5_files_list)
    message('\t  Done! Found ', num_files, ' files')

    message('\t- Starting a parallel compute cluster...\r')
    # Initiate cluster
    cl <- parallel::makeCluster(num_cores)
    doSNOW::registerDoSNOW(cl)
    `%dopar%` <- foreach::`%dopar%`

    # progress bar
    pb <- txtProgressBar(min=1, max=num_files, style=3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)
    message('\t  Done!')

    #loop
    message('\t- Searching for Poly(A) and Poly(T) tails...\r')
    mcoptions <- list(preschedule=TRUE, set.seed=FALSE, cleanup=TRUE)
    data_list <-foreach::foreach(file_path = fast5_files_list,
                                 .combine = 'rbind',
                                 .inorder = FALSE,
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
    message('\t  Done!')

    message('\t- Shutting down the parallel compute cluster...\r')
    close(pb)
    parallel::stopCluster(cl)
    message('\t  Done!')

    message('\t- Making a dataframe of the tail data...\r')
    # Round the tail length
    # df$tail_length <- as.numeric(df$tail_length, digits=2)
    df <- data.frame(data_list, stringsAsFactors = FALSE)
    message('\t  Done!')

    message('\t- Saving the dataframe...\r')
    data.table::fwrite(df, file.path(save_dir, csv_filename, fsep = .Platform$file.sep))

    # unlist the columns of the dataframe
    df <- t(apply(df, 1, unlist))
    # convert to dataframe again
    df <- data.frame(df, stringsAsFactors = FALSE)
    # change data types of colummns from character to numeric
    df$tail_length <- as.numeric(df$tail_length)
    df$tail_start <- as.numeric(df$tail_start)
    df$tail_end <- as.numeric(df$tail_end)
    df$samples_per_nt <- as.numeric(df$samples_per_nt)
    message('\t  Done!')

    return(df)
}
