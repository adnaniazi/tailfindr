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

    # Read read_ids from Fast5s
    message('Step 1 of 3: Finding read IDs from Fast5 files\n')
    df_fast5 <- data.frame(get_fast5_read_ids_parallel_hdf5r(fast5_dir))
    message('\nDone!\n')

    # get read_ids and strand info from Fast5s
    message('Step 2 of 3: Finding forward and reverse strand information from the BAM file\n')
    message('This may take a while; please wait....')
    df_bam <- segregate_forward_and_reverse_strand_reads (alignment_bam_file)
    message('Done!\n')
    df_bam <- dplyr::rename(df_bam, read_id = qname)
    # merge the two dfs
    df <- dplyr::inner_join(df_fast5, df_bam, by='read_id')
    df <- na.omit(df)

    # get to positive and negative strand dfs
    df_polya <- dplyr::filter(df, strand=='+')
    df_polyt <- dplyr::filter(df, strand=='-')


    # read poly(A) tails in forward strands
    # read poly(T) tails in reverse strands


    message('Step 2 of 3: Finding forward and reverse strand information from the BAM file\n')


}
