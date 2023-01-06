library(plyr)
library(readr)
file_vec1 <-c(
  '/export/valenfs/data/processed_data/MinION/2_prime_o_methylation_project/0_rebasecalling/1_AG/x_kmer_dwell_times/kmer_dwell.csv',
  '/export/valenfs/data/processed_data/MinION/2_prime_o_methylation_project/0_rebasecalling/6_cleancap_cap0/x_kmer_dwell_times/kmer_dwell.csv',
  '/export/valenfs/data/processed_data/MinION/2_prime_o_methylation_project/0_rebasecalling/7_zebrafish_IVT_cap0/x_kmer_dwell_times/kmer_dwell.csv',
  '/export/valenfs/data/processed_data/MinION/2_prime_o_methylation_project/0_rebasecalling/8_zebrafish_IVT_cap0/x_kmer_dwell_times/kmer_dwell.csv',
  '/export/valenfs/data/processed_data/MinION/2_prime_o_methylation_project/0_rebasecalling/22_debruijn_ref_rna_teshome/x_kmer_dwell_times/kmer_dwell.csv'
)
file_vec2 <- list.files("/export/valenfs/data/processed_data/MinION/2_prime_o_methylation_project/0_rebasecalling/21_ec_sc_ivt/x_kmer_dwell_times/", pattern="*.csv", full.names=TRUE)
file_vec <- c(file_vec1, file_vec2)

#file_vec <- file_vec[1:3]

dat_csv = ldply(file_vec, read_csv)
dat_csv<-dat_csv[!(dat_csv$cigar=='D'), ]
detach(package:plyr)


# turns out that dataframe entries are too large in number
# to be handled by dplyr. Switching to data table instead

require(data.table) # v1.9.0+
setDT(dat_csv)

dat_csv[, cigar:=NULL] # reomve cigar column
summary <- dat_csv[,lapply(.SD,median),by=ref_kmer]

library(dplyr)
df <- as.data.frame(summary)
names(df) <- c('ref_kmer', 'median')
#dat_csv <- dat_csv %>% group_by(ref_kmer) %>% partition(cluster)
#summary <- dat_csv %>% summarize(median = median(raw_dwell_time, na.rm = TRUE)) %>% collect()


