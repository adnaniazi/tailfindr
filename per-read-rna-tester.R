rm(list=ls())
library(tailfinder)
#folder_path <- '/export/valenfs/data/processed_data/MinION/20180912_max_rna_256_polya/basecalled_data/workspace/pass/0/'
folder_path <- '/export/valenfs/data/processed_data/MinION/20180515_1725_polya_direct_rna_shield/basecalled_data/workspace/pass/30/'
#folder_path <- '/export/valenfs/data/processed_data/MinION/20180515_1725_polya_direct_rna_shield/polya_estimation/gfp_aligned_reads'

file_path <- list.files(path=folder_path, full.names = TRUE)
# a <- find_rna_polya_tail_per_read(file_path[33],
#                                   save_plots=TRUE,
#                                   show_plots=TRUE,
#                                   save_dir='/export/valenfs/data/processed_data/MinION/20180515_1725_polya_direct_rna_shield/',
#                                   plotting_library='rbokeh',
#                                   plot_debug=T)


a <- find_rna_polya_tail_per_read(file_path[61],
                                  save_plots=TRUE,
                                  show_plots=TRUE,
                                  save_dir='/export/valenfs/data/processed_data/MinION/20180515_1725_polya_direct_rna_shield/',
                                  plotting_library='ggplot2',
                                  plot_debug=T)
