#' Find Poly(T) tails in FAST5 reads
#'
#' This function gets a directory of forward-strand FAST5 reads, searches
#' for poly(T) tails in these reads, and then returns a dataframe that contains
#' the data about poly(T) tails found
#'
#' @param fast5_files_list a character string.
#' @param save_dir a character string. Path of the directory in which to save the
#' results csv file and plots (optional)
#' @param csv_file_name a character string. Name of the csv file in which to store the results.
#' @param save_plots A logical [FALSE]. Set it to TRUE if you want to
#' save a plot of poly(T) tail in each read, otherwise set it to FALSE.
#' @param show_plots A logical [FALSE]. Set it to TRUE if you want to display a
#' plot of poly(T) tails in each read, otherwise set it to FALSE. Setting this
#' option to TRUE could cause the algorithm to run very slowly.
#' @param num_cores
#'
#' @return A dataframe containing information about poly(T) tail lengths of each read
#' @export
#'
#' @examples
#' df <- find_cdna_polyt_tails_batch_mclapply('/FORWARD/STRAND/FAST5/FILES/DIRECTORY', 'SAVE/DIR', 'polya-tail-data.csv')
find_cdna_polyt_tails_batch_mclapply <- function(fast5_files_list,
                                                 save_dir,
                                                 csv_file_name,
                                                 save_plots=FALSE,
                                                 show_plots=FALSE,
                                                 num_cores=1){

    num_files <- length(fast5_files_list)
    message('\t- ', num_files, ' forward-strand Fast5 files to search Poly(T)-tails in.')

    # divide data into chunks such that each chunck has 10K reads
    files_per_chunk <- 10000
    num_chunks_floor <- floor(num_files/files_per_chunk)
    num_chunks_ceil <- ceiling(num_files/files_per_chunk)

    message('\t- Divided these ', num_files, ' reads into ', num_chunks_ceil, ' chunks')

    i = 0
    tmp_list <- list()
    while (i < num_chunks_ceil){
        chunk_start <- (i*files_per_chunk) + 1
        chunk_end <- (i+1)*files_per_chunk
        file_indices <- chunk_start:chunk_end
        i <- i + 1
        if (i == num_chunks_ceil){
            chunk_start <- ((i-1)*files_per_chunk) + 1
            chunk_end <- num_files
            file_indices <- chunk_start:chunk_end
        }

        reads_subset = fast5_files_list[file_indices]
        message('\t- Processing chunk ', i, ' of ', num_chunks_ceil)
        tmp <- parallel::mclapply(reads_subset, find_cdna_polyt_tail_per_read,
                                  show_plots=show_plots, save_plots=save_plots,
                                  save_dir=save_dir, mc.cores = num_cores)
        #tmp <- 1
        tmp_list[[i]] <- tmp
        message('\t  Done!')
    }
    # Unlist and make a dataframe
    poly_t_df <- data.frame()
    for(k in tmp_list){
        df <- data.frame(do.call(rbind, k))
        poly_t_df <- rbind(poly_t_df, df)

    }
    poly_t_df$read_id <- as.character(poly_t_df$read_id)
    message('\t- Finished searching all reads for poly(T) tails.')
    return(poly_t_df)
}

