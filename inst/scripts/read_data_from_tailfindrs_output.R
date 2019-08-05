df <- read.csv("/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20190719_DNA_Max_RCA-GFP-2/tailfindr_rca_umi/tails_string.csv",
               stringsAsFactors = FALSE, header = TRUE)

library(dplyr)
df_fq <- dplyr::select(df, fastq_length, content_string) %>%
    dplyr::arrange(fastq_length)

write.table(df_fq,
            file="/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20190719_DNA_Max_RCA-GFP-2/tailfindr_rca_umi/read_length_vs_content_all.txt",
            col.names = FALSE,
            row.names = FALSE
            )


df_fq <- df_fq[df_fq$content_string != "",]
write.table(df_fq,
            file="/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20190719_DNA_Max_RCA-GFP-2/tailfindr_rca_umi/read_length_vs_content_only_complete_cases.txt",
            col.names = FALSE,
            row.names = FALSE
)

