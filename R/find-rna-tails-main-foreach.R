#' Title
#'
#' @param fast5_dir
#' @param save_dir
#' @param poly_a_csv_file_name
#' @param save_plots
#' @param plot_debug
#' @param plotting_library
#' @param num_cores
#'
#' @return
#' @export
#'
#' @examples
find_rna_tails_foreach <- function(fast5_dir,
                                   save_dir,
                                   poly_a_csv_file_name,
                                   save_plots=FALSE,
                                   show_plots=FALSE,
                                   plot_debug=FALSE,
                                   plotting_library='ggplot2',
                                   num_cores=1){

    poly_a_csv_file <- file.path(save_dir, poly_a_csv_file_name)

    message('\t- Searching for all Fast5 files...\r')
    fast5_files_list <- list.files(path=fast5_dir,
                                   pattern="\\.fast5$",
                                   recursive=TRUE,
                                   full.names=TRUE)
    #fast5_files_list <- fast5_files_list[1:2000]
    polya_tails <- find_rna_polya_tails_foreach(fast5_files_list,
                                                save_dir=save_dir,
                                                save_plots=save_plots,
                                                show_plots=show_plots,
                                                plot_debug=plot_debug,
                                                plotting_library=plotting_library,
                                                num_cores=num_cores)

    # include mapping information in the results
    polya_tails <- unique(polya_tails)
    data.table::fwrite(polya_tails, poly_a_csv_file)
    message('Poly(A) tail metadata has been saved in the following file:')
    message(poly_a_csv_file)

    message('Completed step 4 successfully!\n')
}
