#' Find Poly(A) tails in FAST5 reads
#'
#' This function gets a directory of forward-strand FAST5 reads, searches
#' for poly(A) tails in these reads, and then returns a dataframe that contains
#' the data about poly(A) tails found
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
#' df <- find_cdna_polya_tails_batch_mclapply('/FORWARD/STRAND/FAST5/FILES/DIRECTORY', 'SAVE/DIR', 'polya-tail-data.csv')
find_cdna_polya_tails_batch_mclapply <- function(fast5_files_list,
                                                 poly_a_adaptor="GAAGATAGAGCGACAGGCAAGT",
                                                 save_dir,
                                                 csv_file_name,
                                                 save_plots=FALSE,
                                                 show_plots=FALSE,
                                                 num_cores=1){
    num_files <- length(fast5_files_list)
    message('\t- ', num_files, ' forward-strand Fast5 files to search Poly(A)-tails in.')

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
        tmp <- parallel::mclapply(reads_subset, find_cdna_polya_tail_per_read,
                                  poly_a_adaptor=poly_a_adaptor,
                                  show_plots=show_plots, save_plots=save_plots,
                                  save_dir=save_dir, mc.cores = num_cores)
        #tmp <- 1
        tmp_list[[i]] <- tmp
        message('\t  Done!')
    }
    # Unlist and make a dataframe
    poly_a_df <- data.frame()
    for(k in tmp_list){
        df <- data.frame(do.call(rbind, k))
        poly_a_df <- rbind(poly_a_df, df)

    }
    poly_a_df$read_id <- as.character(poly_a_df$read_id)
    message('\t- Finished searching all reads for poly(A) tails.')
    return(poly_a_df)
}

