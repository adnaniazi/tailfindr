library(tailfinder)
fast5_dir <- '~/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20181122_dna_pcr_spikeins/basecalled_data/workspace/pass/0'
alignment_bam_file <- '~/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20181122_dna_pcr_spikeins/alignment_to_gfp/filtered.aln.bam'
save_dir <- '~/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20181122_dna_pcr_spikeins/spikeins_polya_polyt_data'
poly_a_csv_file_name <- 'polya-tails.csv'
poly_t_csv_file_name <- 'polyt-tails.csv'

save_plots <- FALSE
num_cores <- 10
fast5_files_list <- list.files(path=fast5_dir,
                               pattern="\\.fast5$",
                               recursive=TRUE,
                               full.names=TRUE)

fast5_files_list <- fast5_files_list[1:10]

a <- find_cdna_polyt_tails_foreach (fast5_files_list,
                                    poly_a_adaptor="GCGCCGCCCGGACTC",
                                    save_dir,
                                    poly_t_csv_file_name,
                                    save_plots=FALSE,
                                    show_plots=FALSE,
                                    num_cores=1)
