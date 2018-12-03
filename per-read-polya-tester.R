#rm(list=ls())
#file_path <- '/Home/ii/adnann/code/tailfinder/data-raw/cdna-polya-reads/28.fast5'
file_path <- '~/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20181122_dna_pcr_spikeins/basecalled_data/workspace/pass/160/sars_HP_Z240_Tower_Workstation_20181122_FAK23880_MN21607_sequencing_run_all_82961_read_114123_ch_247_strand.fast5'
#folder_path <- '/export/valenfs/data/processed_data/MinION/20180918_max_cdna_256_polya/demultiplexed_gfp_spikeins/poly_a_reads_100_percent_match_threshold/fished_out_reads/_10bp_reads/pass/0/'

#file_path <- list.files(path=folder_path, full.names = TRUE)
a <- find_cdna_polya_tail_per_read(file_path, save_plots=TRUE, show_plots=TRUE, save_dir='~')


