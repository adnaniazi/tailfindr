# rm(list=ls())
# library(tailfinder)
# library(ggplot2)
# devtools::load_all()
#
# file_path <- system.file(package="tailfinder",
#                          "data-raw",
#                          "cdna-polya-reads",
#                          "0.fast5")
#
# read_data <- extract_read_data(file_path)
# norm_data <- z_normalize(read_data$raw_data)
# truncated_data <- truncate_spikes(norm_data, spike_threshold=.1)
#
# df = data.frame(x=c(1:length(read_data$raw_data)),
#                 raw_data=read_data$raw_data,
#                 norm_data=norm_data,
#                 truncated_data=truncated_data,
#                 moves=read_data$moves_sample_wise_vector)
#
# ggplot2::ggplot(data = df, ggplot2::aes(x = x)) +
#     ggplot2::geom_line(ggplot2::aes(y = norm_data), color='red') +
#     ggplot2::geom_line(ggplot2::aes(y = truncated_data), color='blue') +
#     ggplot2::geom_line(ggplot2::aes(y = moves))

rm(list=ls())
#fast5_dir <- '/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20180516_1429_polya_cdna_shield_run1/basecalled_data/workspace/pass/0'
fast5_dir <- '/Users/adnaniazi/Documents/phd/delete/1'
alignment_bam_file <- '/Users/adnaniazi/Documents/phd/delete/bam/sorted.aln.bam'
#alignment_bam_file <- '/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20180516_1429_polya_cdna_shield_run1/alignment_to_genome/aln.bam'

find_cdna_tails(fast5_dir=fast5_dir,
                alignment_bam_file=alignment_bam_file,
                tails='both',
                save_dir='.',
                csv_file_name='.',
                save_plots=FALSE)
