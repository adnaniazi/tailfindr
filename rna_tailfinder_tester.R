library(tailfinder)
fast5_dir <- '/export/valenfs/data/processed_data/MinION/20180515_1725_polya_direct_rna_shield/polya_estimation/gfp_aligned_reads'
save_dir <- '/export/valenfs/data/processed_data/MinION/20180515_1725_polya_direct_rna_shield/polya_estimation'
poly_a_csv_file_name <- 'poly_a_data_adnans_algorithm.csv'
plotting_library <-  'rbokeh'
save_plots <-  TRUE
plot_debug <- TRUE
num_cores <- 120

df2 <- find_rna_tails_foreach(fast5_dir,
                       save_dir,
                       poly_a_csv_file_name,
                       save_plots=save_plots,
                       plot_debug=plot_debug,
                       plotting_library=plotting_library,
                       num_cores=num_cores)


poly_a_csv_file <- file.path(save_dir, poly_a_csv_file_name)
df2 <- read.csv(file=poly_a_csv_file, header=TRUE, sep=",")


df <- read.csv(file='/export/valenfs/data/processed_data/MinION/20180515_1725_polya_direct_rna_shield/polya_estimation/filename_barcodes.csv',
               header=TRUE, sep=",")
df2 <- read.csv(file='/export/valenfs/data/processed_data/MinION/20180515_1725_polya_direct_rna_shield/polya_estimation/poly_a_data_adnans_algorithm.csv',
               header=TRUE, sep=",")


df3 <- dplyr::inner_join(df, df2) %>% na.omit()

df_10 <- dplyr::filter(df3, barcode==10)
df_30 <- dplyr::filter(df3, barcode==30)
df_60 <- dplyr::filter(df3, barcode==60)
df_100 <- dplyr::filter(df3, barcode==100)
df_150 <- dplyr::filter(df3, barcode==150)

p <- ggplot(data=df3, aes(x=poly_a_length_in_nucleotides_1)) +
    stat_density(aes(group = as.factor(barcode), color = as.factor(barcode)),position="identity",geom="line")

p <- ggplot(data=df3, aes(x=as.factor(barcode))) +
    geom_boxplot(aes(y=poly_a_length_in_nucleotides_1)) +
    scale_y_continuous(limits = c(0,200))

p

library(plotly)
ggplotly(p)
