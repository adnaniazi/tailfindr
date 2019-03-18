rm(list=ls())
tails <- data.table::fread('/export/valenfs/data/processed_data/MinION/20180516_1429_polya_cdna_shield_run1/bullshit_delete_later/polya-tails.csv')
png_dir <- '/export/valenfs/data/processed_data/MinION/20180516_1429_polya_cdna_shield_run1/bullshit_delete_later/plots'
plots <- data.frame(list.files(path=png_dir, pattern="\\.png$", recursive=TRUE, full.names=TRUE))
plots <- dplyr::rename(plots, file_path=list.files.path...png_dir..pattern.......png....recursive...TRUE..)
plots$file_path <- as.character(plots$file_path)
tails$file_path <- as.character(tails$file_path)

# reomve png and path before the filename
plots <- dplyr::mutate(plots, file_path=basename(file_path))
plots <- dplyr::mutate(plots, file_path=tools::file_path_sans_ext(file_path))

tails <- dplyr::mutate(tails, file_path=basename(file_path))
tails <- dplyr::mutate(tails, file_path=tools::file_path_sans_ext(file_path))



a <- dplyr::anti_join(tails, plots, by='file_path')
