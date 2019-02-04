#' Finds poly(A) and poly(T) tails in 1D DNA reads
#'
#' This is the main function to find poly(A) and poly(T) lengths in Oxford
#' Nanopore's 1D DNA reads. The function can handle reads that have been
#' basecalled using the standard model or the flip-flop model. Furthermore,
#' it can also handle multifast5 files as well. The function saves a CSV file
#' containing all the tail information, and also returns a tibble containing
#' the same information for further processing by the end-user.
#'
#' @param fast5_dir a character string. Full path of the directory to search the
#' fast5 files in. The direcotry is searched recursively.
#'
#' @param save_dir a character string. Full path of the directory where the CSV
#' file containing the tail lengths should be stored. If save_plots is set to
#' \code{TRUE}, then plots showing the poly(A)/(T) tails are stored within the
#' \code{plots} directory within the \code{save_dir}. This \code{plots}
#' directory is created automatically.
#'
#' @param csv_filename a character string [\code{"tails.csv"}]. Filename of the
#' CSV file in which to store the tail length data
#'
#' @param num_cores a numeric [1]. Num of phyiscal cores to use in processing
#' the data. If you have 4 physical cores in the computer that you are using
#' tailfinder on, then use 3 for \code{num_cores}. Always use 1 less than the
#' number of cores at your disposal.
#'
#' @param save_plots a logical [\code{FALSE}]. If set to \code{TRUE}, a plots
#' directory will be created within the save_dir, and plots showing poly(A) and
#' poly(T) tails on the raw squiggle will be saved in this \code{plots}
#' directory. Creating plots and saving them to the disk is a slow process. So
#' we recommend that you keep this option set to FALSE. If you still want to
#' create plot, we recommend that you run tailfinder on a subset of reads with
#' \code{save_plots} set to \code{TRUE}
#'
#' @param plot_debug_traces a logical [\code{FALSE}]. If set to \code{TRUE},
#' then we will plot debugging information in the plots as well, such as the
#' mean signal, the slope signal, the thresholds, the smoothened signal, etc.
#' We use this option internally to debug our algorithm. This option works only
#' if \code{save_plots} is also set to \code{TRUE}.
#'
#' @param plotting_library a character string [\code{"rbokeh"}]. \code{rbokeh}
#' is the default plotting library that we will use if \code{save_plots} is set
#' to \code{TRUE}. The plots will be saved as \code{<read_id>.html} files in
#' the \code{/save_dir/plots} directory. You can open these HTLM files in any
#' web-browser and interactively view the plots showing the tail region in the
#' raw squiggle. If this option is set to \code{'ggplot2'}, then the polts will
#' be saved as \code{.png} files.
#'
#' @param ... A list. A list of optional parameters. This is currently, reserved
#' for internal use only. By default, DNA reads are assumed to be from a direct
#' cDNA or an amplified cDNA library. However, if the data is from a PCR DNA
#' library, then an additional parameter named \code{data} should be passed with
#' its value set to \code{'pcr-dna'}. This will ensure that the algorithm uses
#' the correct adaptor sequences for the PCR DNA protocol.
#'
#' @return A data tibble is returned containing all the information
#' about the tails found. Always save this returned tibble in a variable (see
#' examples below), otherwise the very long tibble will be printed to the
#' console, which will hang up your rsession.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' library(tailfinder)
#'
#' # 1. Suppose you have 11 cores at your disposal, then you should run tailfinder
#' # on your data as following:
#' df <- find_tails(fast5_dir = '/path/to/fast5/folder/',
#'                  save_dir = '/path/to/a/folder/where/results/are/to/be/stored/',
#'                  csv_filename = 'tails.csv',
#'                  num_cores = 10)
#'
#' # 2. If you also want to save the plots showing the tail on the raw
#' # squiggle using ggplot2 (plots will be save as .png files),
#' # then you should run tailfinder as following:
#' df <- find_tails(fast5_dir = '/path/to/fast5/folder/',
#'                  save_dir = '/path/to/a/folder/where/results/are/to/be/stored/',
#'                  csv_filename = 'tails.csv',
#'                  num_cores = 10,
#'                  save_plots = TRUE,
#'                  plotting_library = 'ggplot2')
#'
#' # 3. If you want to save interactive HTML plots using rbokeh,
#' # then you should run tailfinder as following:
#' df <- find_tails(fast5_dir = '/path/to/fast5/folder/',
#'                  save_dir = '/path/to/a/folder/where/results/are/to/be/stored/',
#'                  csv_filename = 'tails.csv',
#'                  num_cores = 10,
#'                  save_plots = TRUE,
#'                  plotting_library = 'rbokeh')
#'
#' # 4. If you also want to plot debug traces, then you should run tailfinder as
#' # below:
#' df <- find_tails(fast5_dir = '/path/to/fast5/folder/',
#'                  save_dir = '/path/to/a/folder/where/results/are/to/be/stored/',
#'                  csv_filename = 'tails.csv',
#'                  num_cores = 10,
#'                  save_plots = TRUE,
#'                  plot_debug_traces = TRUE,
#'                  plotting_library = 'rbokeh')
#'
#' # N.B.: Making and saving plots is a computationally slow process.
#' # Only generate plots by running tailfinder on a small subset of your reads.
#' }
#'
find_tails <- function(fast5_dir,
                       save_dir,
                       csv_filename = 'tails.csv',
                       num_cores = 1,
                       save_plots = FALSE,
                       plot_debug_traces = FALSE,
                       plotting_library = 'rbokeh',
                       ...) {

    plot_debug <- plot_debug_traces

    # Taking out these parameter from the function parameters list
    # as they may be dangerous for normal users
    show_plots <- FALSE
    if ("data" %in% names(...)) {
        data <- ...$data
    } else {
        data <- 'cdna'
        # data parameter is used only when the experiment_type is dna
    }

    message('\t- Started tailfinder...\r')

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
        message(paste('\t  ', file.path(save_dir, 'plots', fsep = .Platform$file.sep), sep = ''))
    }

    # search for all the fast5 files in the user-specified directory
    message('\t- Searching for all Fast5 files...\r')
    fast5_files_list <- list.files(path = fast5_dir,
                                   pattern = "\\.fast5$",
                                   recursive = TRUE,
                                   full.names = TRUE)
    num_files <- length(fast5_files_list)
    message('\t  Done! Found ', num_files, ' fast5 files.\r')

    # read the first read in the list of reads,
    # and determine all the properties of the data
    message('\t- Analyzing a single read to explore your data...\r')
    type_info <- explore_basecaller_and_fast5type(fast5_files_list[1])
    basecalled_with <- type_info$basecalled_with
    multifast5 <- ifelse(type_info$fast5type == 'multi', TRUE, FALSE)
    experiment_type <- type_info$experiment_type
    read_is_1d <- type_info$read_is_1d
    model <- type_info$model
    if (basecalled_with == 'albacore'){
        message('\t  * The data has been basecalled using Albacore.\r')
    } else {
        message('\t  * The data has been basecalled using Guppy.\r')
    }
    if (model == 'flipflop'){
        message('\t  * Flipflop model was used during basecalling.\r')
    } else {
        message('\t  * Standard model was used during basecalling.\r')
    }
    if (multifast5){
        message('\t  * The reads are packed in multi-fast5 file(s).\r')
    } else {
        message('\t  * Every read is in a single fast5 file of its own.\r')
    }
    if (experiment_type == 'rna'){
        message('\t  * The experiment type is RNA, so we will search\r')
        message('\t    for poly(A) tails.\r')
    } else {
        message('\t  * The experiment type is DNA, so we will search\r')
        message('\t    for both poly(A) and poly(T) tails.\r')
    }
    if (read_is_1d == TRUE){
        message('\t  * The reads are 1D reads.')
    } else {
        message('\t  * The reads are not 1D. Currently, we only support\r')
        message('\t    1D reads. If you believe your reads are 1D, and you are\r')
        message('\t    getting this message erroneously, please feel free\r')
        message('\t    to contact us at adnan.niazi@uib.no. Do not forget to\r')
        message('\t    send us one of the problematic reads so that we can\r')
        message('\t    debug our software, and send you a patch.\r')
        return(0)
    }

    # Make a computer cluster
    message('\t- Starting a parallel compute cluster...\r')
    cl <- parallel::makeCluster(num_cores, outfile = '')
    on.exit(parallel::stopCluster(cl))

    doSNOW::registerDoSNOW(cl)
    `%dopar%` <- foreach::`%dopar%`
    `%do%` <- foreach::`%do%`
    message('\t  Done!')
    mcoptions <- list(preschedule = TRUE, set.seed = FALSE, cleanup = FALSE)

    # If the fast5 are multifast5, then build an index of all the reads within these files
    if (multifast5) {
        message('\t- Discovering reads in the ',  num_files, ' multifast5 files...\r')
        read_id_fast5_file <- dplyr::tibble(read_id = character(),
                                            fast5_file = character())
        for (fast5_file in fast5_files_list) {
            f5_obj <- hdf5r::H5File$new(fast5_file, mode = 'r')
            f5_tree <- f5_obj$ls(recursive = F)
            f5_tree <- f5_tree$name
            f5_tree <- dplyr::mutate(dplyr::tbl_df(f5_tree), fast5_file = fast5_file)
            f5_tree <- dplyr::rename(f5_tree, read_id = value)
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
        if (experiment_type == 'dna') {
            message('\t- Searching for Poly(A) and Poly(T) tails...\r')
        } else {
            message('\t- Searching for Poly(A) tails...\r')
        }
        counter <- 0
        result <- list()
        for(chunk in c(1:total_chunks)){
            # divide data in chunks
            if(chunk == total_chunks) {
                read_id_fast5_file_subset <- read_id_fast5_file[((counter*files_per_chunk)+1):total_files]
                # if the last chunk has only one read then just duplicate this
                # one read so that the progress works and does not throw an error
                if (length(read_id_fast5_file_subset) == 1) {
                    read_id_fast5_file_subset[[2]] <- read_id_fast5_file_subset
                }
            } else {
                read_id_fast5_file_subset <- read_id_fast5_file[((counter*files_per_chunk)+1):((counter+1)*files_per_chunk)]
            }
            counter <- counter + 1
            message(paste('\t  Processing chunk ', chunk, ' of ', total_chunks, '\r', sep = ''))

            # progress bar
            pb <- txtProgressBar(min = 1,
                                 max = length(read_id_fast5_file_subset),
                                 style = 3)
            progress <- function(n) setTxtProgressBar(pb, n)
            opts <- list(progress = progress)

            # foreach loop
            if (experiment_type == 'dna') {
                data_list <- foreach::foreach(riff = read_id_fast5_file_subset,
                                              .combine = 'rbind',
                                              .inorder = FALSE,
                                              .errorhandling = 'pass',
                                              .options.snow = opts,
                                              .options.multicore = mcoptions) %dopar% {
                                                  tryCatch({
                                                      find_dna_tail_per_read(read_id_fast5_file = riff,
                                                                             file_path = NA,
                                                                             data = data,
                                                                             save_plots = save_plots,
                                                                             show_plots = show_plots,
                                                                             plot_debug = plot_debug,
                                                                             save_dir = save_dir,
                                                                             plotting_library = plotting_library,
                                                                             multifast5 = multifast5,
                                                                             basecalled_with = basecalled_with,
                                                                             model = model)
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
            } else {
                data_list <- foreach::foreach(riff = read_id_fast5_file_subset,
                                              .combine = 'rbind',
                                              .inorder = FALSE,
                                              .errorhandling = 'pass',
                                              .options.snow = opts,
                                              .options.multicore = mcoptions) %dopar% {
                                                  tryCatch({
                                                      find_rna_polya_tail_per_read(file_path = NA,
                                                                                   read_id_fast5_file = riff,
                                                                                   multifast5 = multifast5,
                                                                                   basecalled_with = basecalled_with,
                                                                                   save_plots = save_plots,
                                                                                   show_plots = show_plots,
                                                                                   save_dir = save_dir,
                                                                                   plotting_library = plotting_library,
                                                                                   plot_debug = plot_debug)
                                                  },
                                                  error=function(e){
                                                      ls <- list(read_id = NA,
                                                                 tail_start = NA,
                                                                 tail_end = NA,
                                                                 samples_per_nt = NA,
                                                                 tail_length = NA,
                                                                 polya_fastq = NA,
                                                                 file_path = file_path)
                                                  })
                                              }

            }
            result[[chunk]] <- data_list
        }
    } else if (!multifast5) {
        # Split the data into chunks
        files_per_chunk <- 4000
        total_files <- length(fast5_files_list)
        total_chunks <- ceiling(total_files/files_per_chunk)

        counter <- 0
        result <- list()

        if (experiment_type == 'dna') {
            message('\t- Searching for Poly(A) and Poly(T) tails...\r')
        } else {
            message('\t- Searching for Poly(A) tails...\r')
        }

        for(chunk in c(1:total_chunks)) {
            if(chunk == total_chunks) {
                fast5_files_subset <- fast5_files_list[((counter*files_per_chunk)+1):total_files]
                # if the last chunk has only one read then just duplicate this
                # one read so that the progress works and does not throw an error
                if (length(fast5_files_subset) == 1) {
                    fast5_files_subset[[2]] <- fast5_files_subset
                }
            } else {
                fast5_files_subset <- fast5_files_list[((counter*files_per_chunk)+1):((counter+1)*files_per_chunk)]
            }
            counter <- counter + 1

            # progress bar
            message(paste('\t  Processing chunk ', chunk, ' of ', total_chunks, '\r', sep=''))
            pb <- txtProgressBar(min = 1,
                                 max = length(fast5_files_subset),
                                 style = 3)
            progress <- function(n) setTxtProgressBar(pb, n)
            opts <- list(progress = progress)

            # foreach loop
            if (experiment_type == 'dna') {
                data_list <- foreach::foreach(file_path = fast5_files_subset,
                                              .combine = 'rbind',
                                              .inorder = FALSE,
                                              .options.snow = opts,
                                              .options.multicore = mcoptions) %dopar% {
                                                 tryCatch({
                                                     find_dna_tail_per_read(file_path = file_path,
                                                                            data = data,
                                                                            save_plots = save_plots,
                                                                            show_plots = show_plots,
                                                                            plot_debug = plot_debug,
                                                                            save_dir = save_dir,
                                                                            plotting_library = plotting_library,
                                                                            multifast5 = multifast5,
                                                                            basecalled_with = basecalled_with,
                                                                            model = model)
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
            } else {
                data_list <-foreach::foreach(file_path = fast5_files_subset,
                                             .combine = 'rbind',
                                             .inorder = FALSE,
                                             .options.snow = opts,
                                             .options.multicore = mcoptions) %dopar% {
                                                 tryCatch({
                                                     find_rna_polya_tail_per_read(file_path = file_path,
                                                                                  read_id_fast5_file = NA,
                                                                                  multifast5 = multifast5,
                                                                                  basecalled_with = basecalled_with,
                                                                                  save_plots = save_plots,
                                                                                  show_plots = show_plots,
                                                                                  save_dir = save_dir,
                                                                                  plotting_library = plotting_library,
                                                                                  plot_debug = plot_debug)
                                                 },
                                                 error=function(e){
                                                     ls <- list(read_id = NA,
                                                                tail_start = NA,
                                                                tail_end = NA,
                                                                samples_per_nt = NA,
                                                                tail_length = NA,
                                                                polya_fastq = NA,
                                                                file_path = file_path)
                                                 })
                                             }
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
