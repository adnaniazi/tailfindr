library(dplyr)

df <- read.csv(file = "/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20190829_Max_DNA_RCA-3/tailfindr_rca_umi/tails_string_0.8_threshold.csv",
         header = TRUE,
         stringsAsFactors = FALSE)

df <- df %>% dplyr::select(fastq_length, content_string)

df <- df[df$content_string != '', ]
df <- df %>% arrange(fastq_length)

write.table(df,"/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20190829_Max_DNA_RCA-3/tailfindr_rca_umi/tails_string_0.8_threshold.ctxt",
            sep = "\t",
            row.names = FALSE)

