rm(list=ls())
#file_path <- '/Home/ii/adnann/code/tailfinder/data-raw/cdna-polyt-reads/1.fast5'
#file_path <- '/export/valenfs/data/processed_data/MinION/20180516_1429_polya_cdna_shield_run1/basecalled_data/workspace/pass/0/sars_HP_Z240_Tower_Workstation_20180516_FAH88184_MN21607_sequencing_run_shield1_97778_read_120_ch_160_strand.fast5'
#folder_path <- '~/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20181122_dna_pcr_spikeins/basecalled_data/workspace/pass/104/sars_HP_Z240_Tower_Workstation_20181122_FAK23880_MN21607_sequencing_run_all_82961_read_39355_ch_108_strand.fast5'
#folder_path <- '/export/valenfs/data/processed_data/MinION/20180516_1429_polya_cdna_shield_run1/demultiplexed_gfp_spikeins/poly_t_reads_100_percent_match_threshold/fished_out_reads/_10bp_reads/pass/0/'

#file_path <- list.files(path=folder_path, full.names = TRUE)

#file_path <- '/Users/adnaniazi/Documents/phd/code/tailfinder/data-raw/cdna-polyt-reads/52.fast5'
file_path <- '/Users/adnaniazi/mnt/kjempetuja//export/valenfs/data/processed_data/MinION/20181122_dna_pcr_spikeins/polya_estimation/gfp_aligned_reads/poly_t_reads/sars_HP_Z240_Tower_Workstation_20181122_FAK23880_MN21607_sequencing_run_all_82961_read_100037_ch_264_strand.fast5'
a <- find_cdna_polyt_tail_per_read(file_path, save_plots=FALSE, show_plots=TRUE, save_dir='~')

