rm(list=ls())
file_path <- '/Home/ii/adnann/code/tailfinder/data-raw/cdna-polya-reads/3.fast5'
a <- find_cdna_polya_tail_per_read(file_path, save_plots=TRUE, show_plots=TRUE, save_dir='~')
