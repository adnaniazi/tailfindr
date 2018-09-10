#' Title
#'
#' @param fast5_dir
#' @param alignment_bam_file
#' @param tails
#' @param save_dir
#' @param csv_file_name
#' @param save_plots
#'
#' @return
#' @export
#'
#' @examples
find_cdna_tails <- function(fast5_dir, alignment_bam_file, tails='both',
                            save_dir, csv_file_name, save_plots=FALSE){

    ls <- get_fast_read_ids_parallel_hdf5r(fast5_dir)
}
