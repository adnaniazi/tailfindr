get_fast5_read_ids_hdf5r_foreach <- function(fast5_dir, num_cores=1){
    # recursively search for fast5 files
    message('\t- Searching for all Fast5 files...\r')
    fast5_files_list <- list.files(path=fast5_dir,
                                   pattern="\\.fast5$",
                                   recursive=TRUE,
                                   full.names=TRUE)
    message('\t  Done!')

    # Initiate cluster
    message('\t- Starting a parallel cluster...\r')
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
    message('\t- Extracting read IDs...\r')
    mcoptions <- list(preschedule=FALSE, set.seed=FALSE, cleanup=TRUE)
    list_data <- foreach::foreach(file_path = fast5_files_list,
                                  .combine='rbind',
                                  .options.snow=opts,
                                  .options.multicore = mcoptions) %dopar% {
        tryCatch({
            read_id_path <- get_fast5_read_id_hdf5r(file_path)
        },
        error=function(e){
            cat("ERROR :",conditionMessage(e), "\n")
            read_id_path <- list(read_id=NA, file_path=file_path)
        })
    }
    close(pb)
    parallel::stopCluster(cl)

    # cleanup
    df_fast5 <- data.frame(list_data)
    df_fast5$read_id <- as.character(df_fast5$read_id)
    df_fast5$file_path <- as.character(df_fast5$file_path)
    df_fast5 <- na.omit(df_fast5)

    return(df_fast5)
}

