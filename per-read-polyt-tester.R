rm(list=ls())
#file_path <- '/Home/ii/adnann/code/tailfinder/data-raw/cdna-polyt-reads/1.fast5'
file_path <- '/Users/adnaniazi/Documents/phd/code/tailfinder/data-raw/cdna-polyt-reads/53.fast5'
find_cdna_polyt_tail_per_read(file_path, save_plots=FALSE, show_plots=TRUE, save_dir='~')

