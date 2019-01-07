library(tailfinder)

# fast5_dir <- '/export/valenfs/data/processed_data/MinION/20180918_max_cdna_256_polya/gfp_aligned_reads/poly_t_reads/pass/0'
# save_dir <- '/export/valenfs/data/processed_data/MinION/20180918_max_cdna_256_polya/gfp_aligned_reads'
# data <- 'cdna'

fast5_dir <- '/export/valenfs/data/processed_data/MinION/20181122_dna_pcr_spikeins/basecalled_data/workspace/pass/0'
save_dir <- '/export/valenfs/data/processed_data/MinION/20181122_dna_pcr_spikeins/delete_later'
data <- 'cdna'

csv_filename <- 'testing.csv'
save_plots=FALSE
show_plots=FALSE
plot_debug=FALSE
plotting_library='ggplot2'
num_cores <- 120

a <- find_dna_tails(fast5_dir=fast5_dir,
               data=data,
               save_dir=save_dir,
               csv_filename=csv_filename,
               num_cores=num_cores,
               save_plots=save_plots,
               show_plots=show_plots,
               plot_debug=plot_debug,
               plotting_library='rbokeh')


