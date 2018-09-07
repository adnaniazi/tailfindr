extract_read_ids_from_fast5_files_hdf5r <- function(fast5_dir){
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
    ls2<-foreach(file_path = fast5_files_list, .combine='rbind') %dopar% {
        tryCatch({
            read_id_path <- get_fast5_read_id(file_path)
        },
        error=function(e){
            cat("ERROR :",conditionMessage(e), "\n")
            print(e)
        })
    }
    parallel::stopCluster(cl)
    return(ls2)
}



#
# library(hdf5r)
# library(foreach)
# library(doParallel)
# fast5_dir <- '/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20180516_1429_polya_cdna_shield_run1/basecalled_data/workspace/pass/0'
# ls <- extract_read_ids_from_fast5_files_hdf5r(fast5_dir)
