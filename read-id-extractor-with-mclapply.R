fast5_dir <- '/export/valenfs/data/processed_data/MinION/20180516_1429_polya_cdna_shield_run1/basecalled_data/workspace/pass/delete/0'
#a <- get_fast5_read_ids_mclapply_hdf5r(fast5_dir, num_cores=139)


fast5_files_list <- list.files(path=fast5_dir,
                               pattern="\\.fast5$",
                               recursive=TRUE,
                               full.names=TRUE)

b <- find_cdna_polya_tails_batch_mclapply(fast5_files_list,
                                          save_dir='/export/valenfs/data/processed_data/MinION/20180516_1429_polya_cdna_shield_run1/basecalled_data/workspace/pass/delete',
                                          csv_file_name='poly-a-new.csv',
                                          save_plots=FALSE,
                                          show_plots=FALSE,
                                          num_cores=139)

c <- find_cdna_polyt_tails_batch_mclapply(fast5_files_list,
                                          save_dir='/export/valenfs/data/processed_data/MinION/20180516_1429_polya_cdna_shield_run1/basecalled_data/workspace/pass/delete',
                                          csv_file_name='poly-a-new.csv',
                                          save_plots=FALSE,
                                          show_plots=FALSE,
                                          num_cores=139)
