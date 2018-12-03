#' Title
#'
#' @param fast5_dir
#' @param alignment_bam_file
#' @param tails
#' @param save_dir
#' @param poly_a_csv_file_name
#' @param poly_t_csv_file_name
#' @param save_plots
#'
#' @return
#' @export
#'
#' @examples
find_cdna_tails <- function(fast5_dir, alignment_bam_file,
                            tails='both',
                            poly_a_adaptor="GAAGATAGAGCGACAGGCAAGT",
                            poly_t_adaptor="ACTTGCCTGTCGCTCTATCTTC",
                            save_dir,
                            poly_a_csv_file_name,
                            poly_t_csv_file_name,
                            save_plots=FALSE, num_cores=1){

    polya_cadidates_file_path <- file.path(save_dir, 'candidate-poly-a-reads.csv')
    polyt_cadidates_file_path <- file.path(save_dir, 'candidate-poly-t-reads.csv')

    if (!file.exists(polya_cadidates_file_path) &
        !file.exists(polyt_cadidates_file_path)){
        # Read read_ids from Fast5s
        message('Step 1: Extracting read IDs from Fast5 files\r')
        df_fast5 <- get_fast5_read_ids_mclapply_hdf5r(fast5_dir, num_cores=num_cores)
        message('Completed step 1 successfully!\n')

        # get read_ids and strand info from Fast5s
        message('Step 2: Finding forward and reverse strand information from the BAM file\r')
        message('This may take a while depending on the bam file size. Please wait...\r')
        df_bam <- read_bam(alignment_bam_file)

        # merge the two dfs
        df <- dplyr::inner_join(df_fast5, df_bam, by='read_id')
        df <- unique(df)

        # get positive and negative strand data frames
        df_polya <- dplyr::filter(df, strand=='+') # forward strand contains poly(A) tails in cDNA
        df_polyt <- dplyr::filter(df, strand=='-') # reverse strand contains poly(T) tails in cDNA


        # write to file
        data.table::fwrite(df_polya, polya_cadidates_file_path)
        data.table::fwrite(df_polyt, polyt_cadidates_file_path)
        message('Completed step 2 successfully!\n')
    }else {
        # load files from disk
        message('Step 1: Extracting read IDs from Fast5 files (SKIPPED)\r')
        message('Step 2: Finding forward and reverse strand information from the BAM file (SKIPPED)\n')
        message('Why did we skip step 1 and step 2?')
        message('Step 1 and 2 take a lot of time to execute. We have already stored the results')
        message('of these two steps in the following two files:\n')

        message(polya_cadidates_file_path)
        message(polyt_cadidates_file_path)

        message('\nBecause we found these files present in your save_dir, therefore, we are using')
        message('them instead of computing them all over again. That is why we have skipped Step 1 and 2.\n')

        message('If your Fast5 files, or the BAM file has changed in anyway, then you must manually\r')
        message('delete the above two files and re-run the find_cdna_tails() function again.\n')

        df_polya <- data.table::fread(polya_cadidates_file_path)
        df_polyt <- data.table::fread(polyt_cadidates_file_path)
    }

    # find poly(A) and poly(T) tails
    if (tails == 'polyA' | tails == 'both'){
        message('Step 3: Finding Poly(A) tails in forward strand reads\r')
        polya_tails <- find_cdna_polya_tails_batch_mclapply(df_polya$file_path,
                                                            poly_a_adaptor=poly_a_adaptor,
                                                            save_dir=save_dir,
                                                            csv_file_name=csv_file_name,
                                                            save_plots=save_plots,
                                                            show_plots=FALSE,
                                                            num_cores=num_cores)
        # include mapping information in the results
        polya_tails <- dplyr::inner_join(polya_tails,
                                         dplyr::select(df_polya, bam_mapping_quality, read_id),
                                         by='read_id')
        polya_tails <- unique(polya_tails)
        data.table::fwrite(polya_tails, file.path(save_dir, poly_a_csv_file_name))
        message('Poly(A) tail metadata has been saved in the following file:')
        message(file.path(save_dir, poly_a_csv_file_name))

        message('Completed step 3 successfully!\n')

    }

    if (tails == 'polyT' | tails == 'both'){
        if (tails == 'both'){
            message('Step 4: Finding Poly(T) tails in reverse strand reads')
        } else {
            message('Step 3: Finding Poly(T) tails in reverse strand reads')
        }
        polyt_tails <- find_cdna_polyt_tails_batch_mclapply(df_polyt$file_path,
                                                            poly_t_adaptor=poly_t_adaptor,
                                                            save_dir=save_dir,
                                                            csv_file_name=csv_file_name,
                                                            save_plots=save_plots,
                                                            show_plots=FALSE,
                                                            num_cores=num_cores)
        # include mapping information in the results
        polyt_tails <- dplyr::inner_join(polyt_tails,
                                         dplyr::select(df_polyt, bam_mapping_quality, read_id),
                                         by='read_id')
        polyt_tails <- unique(polyt_tails)
        data.table::fwrite(polyt_tails, file.path(save_dir, poly_t_csv_file_name))
        message('Poly(T) tail metadata has been saved in the following file:')
        message(file.path(save_dir, poly_t_csv_file_name))

        if (tails == 'both'){
            message('Completed step 4 successfully!\n')
        } else {
            message('Completed step 3 successfully!\n')
        }
    }

}
