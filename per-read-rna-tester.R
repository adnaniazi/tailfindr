rm(list=ls())
library(tailfinder)
#folder_path <- '/export/valenfs/data/processed_data/MinION/20180912_max_rna_256_polya/basecalled_data/workspace/pass/0/'
folder_path <- '~/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20180515_1725_polya_direct_rna_shield/basecalled_data/workspace/pass/30/'
#folder_path <- '/export/valenfs/data/processed_data/MinION/20180515_1725_polya_direct_rna_shield/polya_estimation/gfp_aligned_reads'

#file_path <- list.files(path=folder_path, full.names = TRUE)
# a <- find_rna_polya_tail_per_read(file_path[33],
#                                   save_plots=TRUE,
#                                   show_plots=TRUE,
#                                   save_dir='/export/valenfs/data/processed_data/MinION/20180515_1725_polya_direct_rna_shield/',
#                                   plotting_library='rbokeh',
#                                   plot_debug=T)


a <- find_rna_polya_tail_per_read('/export/valenfs/data/processed_data/MinION/sars_HP_Z240_Tower_Workstation_20180912_FAJ15086_MN21607_sequencing_run_1209_92041_read_4_ch_65_strand.fast5',
                                  save_plots=FALSE,
                                  show_plots=TRUE,
                                  save_dir='~/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20180515_1725_polya_direct_rna_shield',
                                  plotting_library='ggplot2',
                                  plot_debug=T)
