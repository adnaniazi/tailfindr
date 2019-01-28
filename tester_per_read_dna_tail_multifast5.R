rm(list = ls())
library(hdf5r)
library(tailfinder)

fast5_file <- '/export/valenfs/data/processed_data/MinION/20181219_max_dna_pcrspikes/renanalysis_with_guppy_flipflop/few_4000_per_fast5_reads/1/batch_0.fast5'
read_id_fast5_file <- dplyr::tibble(read_id=character(), fast5_file=character())
for (fast5_file in fast5_file) {
    f5_obj <- hdf5r::H5File$new(fast5_file, mode='r')
    f5_tree <- f5_obj$ls(recursive=F)
    f5_tree <- f5_tree$name
    f5_tree <- dplyr::mutate(dplyr::tbl_df(f5_tree), fast5_file=fast5_file)
    f5_tree <- dplyr::rename(f5_tree, read_id=value)
    read_id_fast5_file <- rbind(read_id_fast5_file, f5_tree)
    f5_obj$close_all()
}
message('\t  Done! Found ', nrow(read_id_fast5_file), ' reads\r')
# convert the data frame to list with rows as elements of the list
read_id_fast5_file <- split(read_id_fast5_file, seq(nrow(read_id_fast5_file)))


data='pcr-dna'
save_dir <- '/export/valenfs/data/processed_data/MinION/20181219_max_dna_pcrspikes/renanalysis_with_guppy_flipflop'
num_cores <- 10
save_plots=F
show_plots=F
plot_debug=F
multifast5=T
basecalled_with_flipflop=T
plotting_library='rbokeh'
read_index = 1


a <- find_dna_tail_per_read(read_id_fast5_file=read_id_fast5_file[[15]],
                            file_path=NA,
                            data=data,
                            save_plots=save_plots,
                            show_plots=show_plots,
                            plot_debug=plot_debug,
                            save_dir=save_dir,
                            plotting_library=plotting_library,
                            multifast5,
                            basecalled_with_flipflop)
event_length_vector <- sort(event_length_vector)
event_length_vector <- event_length_vector[1:ceiling(length(event_length_vector)*0.95)]
hist(event_length_vector, breaks = seq(from=0, to=max(event_length_vector), by=2))

