rm(list = ls())
library(tailfinder)
files <- list.files('/export/valenfs/data/processed_data/MinION/20181219_max_dna_pcrspikes/renanalysis_with_guppy_flipflop/basecalled_data/workspace/0',
                    recursive = T,
                    full.names = T,
                    pattern = '*.fast5')
# df_ff <- data.frame(files)
# which(grepl('sars_HP_Z240_Tower_Workstation_20181219_FAK36778_MN21607_sequencing_run_PCR_48535_read_64665_ch_69_strand.fast5', df_ff$files))
# extract_read_data_hdf5r_flipflop(files[13], plot_debug=T)
#

path_ff <- '/export/valenfs/data/processed_data/MinION/20181219_max_dna_pcrspikes/renanalysis_with_guppy_flipflop/basecalled_data/workspace//54/sars_HP_Z240_Tower_Workstation_20181219_FAK36778_MN21607_sequencing_run_PCR_48535_read_64665_ch_69_strand.fast5'
path_al <- '/export/valenfs/data/processed_data/MinION/20181219_max_dna_pcrspikes/basecalled_data/workspace/pass/36/sars_HP_Z240_Tower_Workstation_20181219_FAK36778_MN21607_sequencing_run_PCR_48535_read_64665_ch_69_strand.fast5'
data_ff <- extract_read_data_hdf5r_flipflop(path_ff, plot_debug=T)
data_al <- extract_read_data_hdf5r(path_al, plot_debug=T)

raw_data <- (data_ff$raw_data - mean(data_ff$raw_data))/ sd(data_ff$raw_data)
library(rbokeh)
df <- data.frame(x=c(1:length(data_ff$raw_data)),
                 raw_data=raw_data,
                 moves_ff=data_ff$moves_sample_wise_vector-5,
                 moves_al=data_al$moves_sample_wise_vector-6)


p <- rbokeh::figure(data=df, width=1000, height=600, legend_location="top_left")
p <- rbokeh::ly_lines(p, raw_data, width=1.5, color='#4040a1', legend = "Raw data")
p <- rbokeh::ly_lines(p, moves_ff, color='green', width=1, legend = "Moves Flip Flop")
p <- rbokeh::ly_lines(p, moves_al, color='orange', width=1, legend = "Moves Albacore")
p

