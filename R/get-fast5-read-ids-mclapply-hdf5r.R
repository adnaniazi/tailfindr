#' Title
#'
#' @param fast5_dir
#' @param num_cores
#'
#' @return
#' @export
#'
#' @examples
get_fast5_read_ids_mclapply_hdf5r <- function(fast5_dir, num_cores=1){
    # recursively search for fast5 files
    message('\t- Searching for all Fast5 files...\r')
    fast5_files_list <- list.files(path=fast5_dir,
                                   pattern="\\.fast5$",
                                   recursive=TRUE,
                                   full.names=TRUE)
    num_files <- length(fast5_files_list)
    message('\t  Found ', num_files, ' Fast5 files.')

    # divide data into chunks such that each chunck has 10K reads
    files_per_chunk <- 10000
    num_chunks_floor <- floor(num_files/files_per_chunk)
    num_chunks_ceil <- ceiling(num_files/files_per_chunk)

    message('\t- Divided ', num_files, ' reads into ', num_chunks_ceil, ' chunks')

    i = 0
    read_data <- list()
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
        message('\t  Processing chunk ', i, ' of ', num_chunks_ceil)
        tmp <- parallel::mclapply(reads_subset, get_fast5_read_id_hdf5r, mc.cores = num_cores)
        #tmp <- 1
        read_data[[i]] <- tmp

    }
    # Unlist and make a dataframe
    read_data_df <- data.frame()
    for(k in read_data){
        df <- data.frame(do.call(rbind, k))
        read_data_df <- rbind(read_data_df, df)

    }
    read_data_df$read_id <- as.character(read_data_df$read_id)
    return(unique(read_data_df))
}
