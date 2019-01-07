library(tailfinder)
polya_folder <- '/export/valenfs/data/processed_data/MinION/20180918_max_cdna_256_polya/gfp_aligned_reads/poly_a_reads'
polyt_folder <- '/export/valenfs/data/processed_data/MinION/20180918_max_cdna_256_polya/gfp_aligned_reads/poly_t_reads/pass/0'
fast5_files_list_polya <- list.files(path=polya_folder,
                                     pattern="\\.fast5$",
                                     recursive=TRUE,
                                     full.names=TRUE)
fast5_files_list_polyt <- list.files(path=polyt_folder,
                                     pattern="\\.fast5$",
                                     recursive=TRUE,
                                     full.names=TRUE)
save_plots=FALSE
show_plots=TRUE
plot_debug=TRUE
save_dir='~'
plotting_library='rbokeh'

find_dna_tail_per_read(fast5_files_list_polya[512],
                       save_plots=save_plots,
                       show_plots=show_plots,
                       plot_debug=plot_debug,
                       save_dir=save_dir,
                       plotting_library=plotting_library)

