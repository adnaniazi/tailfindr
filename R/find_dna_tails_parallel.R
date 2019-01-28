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
#' @param multifast5
#' @param basecalled_with_flipflop
#'
#' @return
#' @export
#'
#' @examples
find_dna_tails_parallel <- function(fast5_dir,
                           data='cdna',
                           save_dir=NA,
                           csv_filename='dna_tails.csv',
                           num_cores=1,
                           save_plots=F,
                           show_plots=F,
                           plot_debug=F,
                           multifast5=F,
                           basecalled_with_flipflop=F,
                           plotting_library='rbokeh'){

    if (multifast5 & !basecalled_with_flipflop) {
        message('Currently, we do not support multifast5 files that have been basecalled with the Albacore algorithm.\r')
        message('To use this tool, convert your files into the legacy format, i.e., one read per fast5 file.\r')
        message('You can use the multi_to_single_fast5 script in the ont_fast5_api to convert multifast5 read into single reads:\r')
        message('https://github.com/nanoporetech/ont_fast5_api/blob/master/README.rst\r')
        message('Quitting for now!\r')
        return(0)
    }

    # Try to create the save directory
    if (!dir.exists(file.path(save_dir))) {
        message('\t- Save dir does not exist. Trying to create it...\r')
        tryCatch({
            dir.create(file.path(save_dir, fsep = .Platform$file.sep))
            message('\t  Done!')
        },
        error=function(e){
            message('\t  Failed to create the save dir. Results will be stored in the "~/" directory instead\r')
            save_dir <- '~/'
        })
    }

    # Create a sub-direcotry to save all the plots
    if (save_plots){
        message('\t- Creating a sub-directory to save the plots in...\r')
        dir.create(file.path(save_dir, 'plots', fsep = .Platform$file.sep))
        message('\t  Done! All plots will be saved in the following direcotry:\r')
        message(paste('\t  ', file.path(save_dir, 'plots', fsep = .Platform$file.sep), sep=''))
    }

    # search for all the fast5 files in the user-specified directory
    message('\t- Searching for all Fast5 files...\r')
    fast5_files_list <- list.files(path=fast5_dir,
                                   pattern="\\.fast5$",
                                   recursive=TRUE,
                                   full.names=TRUE)
    num_files <- length(fast5_files_list)
    message('\t  Done! Found ', num_files, ' fast5 files.\r')

    # Initiate cluster
    message('\t- Starting a parallel compute cluster...\r')
    cl <- parallel::makeCluster(num_cores, outfile='')
    on.exit(parallel::stopCluster(cl))

    doSNOW::registerDoSNOW(cl)
    `%dopar%` <- foreach::`%dopar%`
    `%do%` <- foreach::`%do%`
    message('\t  Done!')

    # If the fast5 are multifast5, then build an index of all the reads within these files
    if (multifast5 & basecalled_with_flipflop) {
        message('\t- Discovering reads in the ',  num_files, ' multifast5 files found...\r')

        read_id_fast5_file <- dplyr::tibble(read_id=character(), fast5_file=character())
        for (fast5_file in fast5_files_list) {
            f5_obj <- hdf5r::H5File$new(fast5_file, mode='r')
            f5_tree <- f5_obj$ls(recursive=F)
            f5_tree <- f5_tree$name
            f5_tree <- dplyr::mutate(dplyr::tbl_df(f5_tree), fast5_file=fast5_file)
            f5_tree <- dplyr::rename(f5_tree, read_id=value)
            read_id_fast5_file <- rbind(read_id_fast5_file, f5_tree)
            f5_obj$close_all()
        }
        message('\t  Done! Found ', nrow(read_id_fast5_file), ' reads\r')
        # convert the data frame to list with rows as elements of the list
        read_id_fast5_file <- split(read_id_fast5_file, seq(nrow(read_id_fast5_file)))

        # Split the data into chunks
        files_per_chunk <- 4000
        total_files <- length(read_id_fast5_file)
        total_chunks <- ceiling(total_files/files_per_chunk)

        #loop
        message('\t- Searching for Poly(A) and Poly(T) tails...\r')
        counter <- 0
        result <- list()
        for(chunk in c(1:total_chunks)){

            # divide data in chunks
            if(chunk == total_chunks)
                read_id_fast5_file_subset <- read_id_fast5_file[((counter*files_per_chunk)+1):total_files]
            else
                read_id_fast5_file_subset <- read_id_fast5_file[((counter*files_per_chunk)+1):((counter+1)*files_per_chunk)]
            counter <- counter + 1
            message(paste('\t  Processing chunk ', chunk, ' of ', total_chunks, '\r', sep=''))

            # progress bar
            message('\r')
            pb <- txtProgressBar(min=1, max=length(read_id_fast5_file_subset), style=3)
            progress <- function(n) setTxtProgressBar(pb, n)
            opts <- list(progress=progress)
            message('\t  Done!')
            # foreach loop
            mcoptions <- list(preschedule=TRUE, set.seed=FALSE, cleanup=TRUE)
            data_list <- foreach::foreach(riff=read_id_fast5_file_subset,
                                          .combine='rbind',
                                          .inorder=FALSE,
                                          .errorhandling='pass',
                                          .options.snow=opts,
                                          .options.multicore=mcoptions) %dopar% {
                                              tryCatch({
                                                  find_dna_tail_per_read(read_id_fast5_file=riff,
                                                                         file_path=NA,
                                                                         data=data,
                                                                         save_plots=save_plots,
                                                                         show_plots=show_plots,
                                                                         plot_debug=plot_debug,
                                                                         save_dir=save_dir,
                                                                         plotting_library=plotting_library,
                                                                         multifast5,
                                                                         basecalled_with_flipflop)
                                              },
                                              error=function(e){
                                                  ls <- list(read_id = riff$read_id,
                                                             read_type = NA,
                                                             tail_is_valid = NA,
                                                             tail_start = NA,
                                                             tail_end = NA,
                                                             samples_per_nt = NA,
                                                             tail_length = NA,
                                                             file_path = riff$fast5_file,
                                                             has_precise_boundary = NA)
                                              })
                                          }
            result[[chunk]] <- data_list
        }

    } else if (!multifast5 & !basecalled_with_flipflop) {
        # Split the data into chunks
        files_per_chunk <- 4000
        total_files <- length(fast5_files_list)
        total_chunks <- ceiling(total_files/files_per_chunk)

        counter <- 0
        result <- list()
        #loop
        message('\t- Searching for Poly(A) and Poly(T) tails...\r')
        for(chunk in c(1:total_chunks)){
            if(chunk == total_chunks)
                fast5_files_subset <- fast5_files_list[((counter*files_per_chunk)+1):total_files]
            else
                fast5_files_subset <- fast5_files_list[((counter*files_per_chunk)+1):((counter+1)*files_per_chunk)]
            counter <- counter + 1

            # progress bar
            message('\r')

            message(paste('\t  Processing chunk ', chunk, ' of ', total_chunks, '\r', sep=''))
            pb <- txtProgressBar(min=1, max=length(fast5_files_subset), style=3)
            progress <- function(n) setTxtProgressBar(pb, n)
            opts <- list(progress=progress)

            mcoptions <- list(preschedule=TRUE, set.seed=FALSE, cleanup=TRUE)
            data_list <-foreach::foreach(file_path = fast5_files_list,
                                         .combine = 'rbind',
                                         .inorder = FALSE,
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
            result[[chunk]] <- data_list
        }

    }

    message('\t- Formatting the tail data...\r')
    result <- purrr::map(result, function(.x) tibble::as_tibble(.x))
    result <- dplyr::bind_rows(result, .id = "chunk")
    result <- dplyr::select(result, -chunk)
    message('\t  Done!')

    message('\t- Saving the data in the CSV file...\r')
    data.table::fwrite(result, file.path(save_dir, csv_filename, fsep = .Platform$file.sep))
    message('\t  Done!')
    return(result)

}
