rm(list=ls())
file_path <- '/Home/ii/adnann/code/tailfinder/data-raw/cdna-polyt-reads/0.fast5'
a <- find_cdna_polyt_tail_per_read(file_path, save_plots=FALSE, show_plots=TRUE, save_dir='~')
