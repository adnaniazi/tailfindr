rm(list = ls())
library(tailfindr)

df <- read.csv("/Users/adnaniazi/Documents/phd/code/tailfindr/inst/scripts/wholeGFP-hit.csv",
               header = TRUE, stringsAsFactors = FALSE)
df$filepath <-  gsub(pattern = "/Users/max/R/RCA/",
                     replacement = "/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20190506_DNA_Max_RCA-GFP/basecalled_data_guppy_301_single_fast5/",
                     x = df$filepath)


# choose a file of interest

# polyA and polyT mixture with truncated gfp
#filepath <- "/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20190506_DNA_Max_RCA-GFP/basecalled_data_guppy_301_single_fast5/4/b86a1494-ae62-4092-a104-5df8f64c22a1.fast5"

# invalid read
#filepath <- "/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20190506_DNA_Max_RCA-GFP/basecalled_data_guppy_301_single_fast5/13/9e10a4a2-f884-4861-ab6e-b0d73667d93a.fast5"

#polyT (good) with poly(A) at the end
#filepath <- "/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20190506_DNA_Max_RCA-GFP/basecalled_data_guppy_301_single_fast5/14/95c012fd-ed98-41ca-8fe4-2101b16a45f8.fast5"


# polyT (good) with no poly(A) at the end
#filepath <- "/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20190506_DNA_Max_RCA-GFP/basecalled_data_guppy_301_single_fast5/12/50a3c9b6-ba0c-4b37-8eae-30e5b90a6fb4.fast5"

# polyA (bad) gfp not found rc_ologiT not found
#filepath <- "/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20190506_DNA_Max_RCA-GFP/basecalled_data_guppy_301_single_fast5/8/eb7f9ffd-c994-4fc5-b5b3-35ee76cf5228.fast5"

# polya 100nt, polyt10nt, GFP found, rc_gFP not found, both oligoT and rc_oligo T found
#filepath <- "/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20190506_DNA_Max_RCA-GFP/basecalled_data_guppy_301_single_fast5/4/721825bf-8b79-4bfd-ab75-250663317cff.fast5"

# few oligo_T
#filepath <- "/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20190506_DNA_Max_RCA-GFP/basecalled_data_guppy_301_single_fast5/13/a4ffe495-9bf3-4611-b26d-806d7d283fc1.fast5"

#filepath <- "/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20190506_DNA_Max_RCA-GFP/basecalled_data_guppy_301_single_fast5/2/9aea82a6-e9ef-4ca4-9a6f-0d7299731f13.fast5"
#filepath <- "/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20190506_DNA_Max_RCA-GFP/basecalled_data_guppy_301_single_fast5/15/c7f631d9-bdfe-40ab-976b-d2780eb2b63b.fast5"

# good polyA polyT
#filepath <- "/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20190506_DNA_Max_RCA-GFP/basecalled_data_guppy_301_single_fast5/1/f802fbac-a7ec-440b-a84f-12f01e9f15ad.fast5"

# good polyT 20 nt
#filepath <- "/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20190506_DNA_Max_RCA-GFP/basecalled_data_guppy_301_single_fast5/8/3fc7b899-ee97-409f-9b42-e90341721528.fast5"

#filepath <- "/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20190506_DNA_Max_RCA-GFP/basecalled_data_guppy_301_single_fast5/6/c79ab9a9-2e43-413e-a25a-bca2d46151b3.fast5"
#filepath <- "/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20190506_DNA_Max_RCA-GFP/basecalled_data_guppy_301_single_fast5/8/de6e8ad7-1e74-4949-9de7-49a8303071a0.fast5"
#filepath <- "/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20190506_DNA_Max_RCA-GFP/basecalled_data_guppy_301_single_fast5/14/a30da377-5433-49cb-afc7-87d40569b270.fast5"
#filepath <- "/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20190506_DNA_Max_RCA-GFP/basecalled_data_guppy_301_single_fast5/4/a4e068d0-a0b3-4d65-b311-2456d6905e1b.fast5"

# fatal reads
filepath <- "/Users/adnaniazi/Documents/phd/delete_later/rca/0d7eb1fa-1758-4280-8894-78a938e45ddb.fast5"

# Extract read data
read_data <- extract_read_data(file_path = filepath,
                               read_id_fast5_file = NA,
                               plot_debug = FALSE,
                               basecalled_with = 'Guppy',
                               basecall_group = 'Basecall_1D_000',
                               multifast5 = FALSE,
                               model = 'flipflop',
                               plotting_library = 'rbokeh')

fastq <- read_data$fastq
fastq_biostring <- Biostrings::DNAString(fastq)

# define parameters for pairwise alignment
match <- 1
mismatch <- -1
submat <- Biostrings::nucleotideSubstitutionMatrix(match = match,
                                                   mismatch = mismatch)

ot <- Biostrings::DNAString("GCCTGTCGCTCTATCTTC")
rc_ot <- Biostrings::DNAString("GAAGATAGAGCGACAGGC")

# STEP1: Find ot and rc_ot oligos
ot_rc_ot_df <- rca_search_oligo(fastq, ot, rc_ot, submat)

# STEP2: Find poly(A) and poly(T) segments within the read
polyat_df <- rca_locate_polya_polyt_segments(ot_rc_ot_df, fastq_biostring)

# STEP3: Cluster reads
data <- rca_cluster_read_segments(polyat_df)

# STEP4: Create consensus of each cluster
rca_data <- rca_create_consensus_sequence(data)

# STEP 5: find tail lengths
rca_data_2 <- rca_find_tails_per_read(rca_data = rca_data,
                                    read_data = read_data,
                                    file_path = filepath)




# run the final algorithm
library(tailfindr)
df_dna <- find_tails(fast5_dir = "/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20190506_DNA_Max_RCA-GFP/basecalled_data_guppy_301_single_fast5/",
                     save_dir = '~/Downloads/rca',
                     basecall_group = 'Basecall_1D_000',
                     csv_filename = 'tails.csv',
                     num_cores = 4,
                     save_plots = FALSE,
                     dna_datatype = 'rca-cdna')


library(tailfindr)
df_dna <- find_tails(fast5_dir = "/export/valenfs/data/processed_data/MinION/20190506_DNA_Max_RCA-GFP/basecalled_data_guppy_301_single_fast5/",
                     save_dir = '/export/valenfs/data/processed_data/MinION/20190506_DNA_Max_RCA-GFP/tailfindr_rca',
                     basecall_group = 'Basecall_1D_000',
                     csv_filename = 'tails.csv',
                     num_cores = 10,
                     save_plots = FALSE,
                     dna_datatype = 'rca-cdna')

df_dna <- find_tails(fast5_dir = "/Users/adnaniazi/Documents/phd/delete_later/rca3",
                     save_dir = '~/Downloads/rca',
                     basecall_group = 'Basecall_1D_000',
                     csv_filename = 'tails.csv',
                     num_cores = 2,
                     save_plots = FALSE,
                     dna_datatype = 'rca-cdna')


df <- rca_find_tails_main(file_path = "/Users/adnaniazi/Documents/phd/delete_later/rca3/50a3c9b6-ba0c-4b37-8eae-30e5b90a6fb4_2.fast5",
                                read_id_fast5_file = NA,
                                basecall_group = 'Basecall_1D_000',
                                save_plots = FALSE,
                                show_plots = FALSE,
                                plot_debug = FALSE,
                                save_dir = '~/Downloads/rca',
                                plotting_library = 'rbokeh',
                                multifast5 = FALSE,
                                basecalled_with = 'guppy',
                                model = 'flipflop')


# Analysis
df <- read.csv("~/Downloads/rca/tails.csv", header = TRUE, stringsAsFactors = FALSE)
