get_fast_read_ids_parallel_hdf5r <- function(fast5_dir){
    fast5_files_list <- list.files(path = fast5_dir,
                                   pattern = "\\.fast5$",
                                   recursive = TRUE,
                                   full.names = TRUE)

    # detect the number of available cores, and use all but one
    no_cores <- parallel::detectCores() - 1

    # Initiate cluster
    cl <- parallel::makeCluster(no_cores)
    doParallel::registerDoParallel(cl)

    `%dopar%` <- foreach::`%dopar%`

    #loop
    ls2<-foreach::foreach(file_path = fast5_files_list, .combine='rbind') %dopar% {
        tryCatch({
            read_id_path <- get_fast5_read_id(file_path)
        },
        error=function(e){
            cat("ERROR :",conditionMessage(e), "\n")
            read_id_path <- list(read_id=NA, file_path=file_path)
        })
    }
    parallel::stopCluster(cl)
    return(ls2)
}

