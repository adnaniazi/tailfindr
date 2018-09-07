#' Find Poly(A) tails in FAST5 reads
#'
#' This function gets a directory of forward-strand FAST5 reads, searches
#' for poly(A) tails in these reads, and then returns a dataframe that contains
#' the data about poly-(A) tails found
#'
#' @param cdna_polya_fast5_dir a character string. Path of the directory containing forward-strand
#' cDNA reads in FAST5 format. This directory will be searched recursively to
#' discover all FAST5 files
#' @param save_dir a character string. Path of the directory in which to save the
#' results csv file and plots (optional)
#' @param csv_file_name a character string. Name of the csv file in which to store the results.
#' @param save_plots A logical [FALSE]. Set it to TRUE if you want to
#' save a plot of poly(A) tail in each read, otherwise set it to FALSE.
#' @param show_plots A logical [FALSE]. Set it to TRUE if you want to display a
#' plot of poly(A) tails in each read, otherwise set it to FALSE. Setting this
#' option to TRUE could cause the algorithm to run very slowly.
#'
#' @return A dataframe containing information about poly(A) tail lengths of each read
#' @export
#'
#' @examples
#' df <- find_cdna_polya_tails('/FORWARD/STRAND/FAST5/FILES/DIRECTORY', 'SAVE/DIR', 'polya-tail-data.csv')
find_cdna_polya_tails <- function(cdna_polya_fast5_dir,
                                  save_dir,
                                  csv_file_name,
                                  save_plots=FALSE,
                                  show_plots=FALSE){

    # recursively search for fast5 files and make a list of them
    fast5_files_list <- list.files(path = cdna_polya_fast5_dir,
                                   pattern = "\\.fast5$",
                                   recursive = TRUE,
                                   full.names = TRUE)


    # detect the number of available cores, and use all but one
    no_cores <- parallel::detectCores() - 1

    # Initiate cluster
    cl <- parallel::makeCluster(no_cores)
    doParallel::registerDoParallel(cl)

    `%dopar%` <- foreach::`%dopar%`
    `%do%` <- foreach::`%do%`

    #loop
    ls2<-foreach(file_path = fast5_files_list, .combine='rbind') %dopar% {
        tryCatch({
            find_read_tail_cdna_polya(file_path,
                                      show_plots=show_plots,
                                      save_plots=save_plots,
                                      save_dir=save_dir)
        },
        error=function(e){
            cat("ERROR :",conditionMessage(e), "\n")
            print(e)
        })

    }
    parallel::stopCluster(cl)


    return(ls2)

}
