#improved RCA detection script

# load the libraries
rm(list = ls())
library(dplyr)
library(tailfindr)
library(Biostrings)
library(msa)
library(DECIPHER)
#load the csv file
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
filepath <- "/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20190506_DNA_Max_RCA-GFP/basecalled_data_guppy_301_single_fast5/12/50a3c9b6-ba0c-4b37-8eae-30e5b90a6fb4.fast5"

# polyA (bad) gfp not found rc_ologiT not found
#filepath <- "/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20190506_DNA_Max_RCA-GFP/basecalled_data_guppy_301_single_fast5/8/eb7f9ffd-c994-4fc5-b5b3-35ee76cf5228.fast5"

# polya 100nt, polyt10nt, GFP found, rc_gFP not found, both oligoT and rc_oligo T found
#filepath <- "/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20190506_DNA_Max_RCA-GFP/basecalled_data_guppy_301_single_fast5/4/721825bf-8b79-4bfd-ab75-250663317cff.fast5"

# few oligo_T
#filepath <- "/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20190506_DNA_Max_RCA-GFP/basecalled_data_guppy_301_single_fast5/13/a4ffe495-9bf3-4611-b26d-806d7d283fc1.fast5"

#filepath <- "/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20190506_DNA_Max_RCA-GFP/basecalled_data_guppy_301_single_fast5/2/9aea82a6-e9ef-4ca4-9a6f-0d7299731f13.fast5"
#filepath <- "/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20190506_DNA_Max_RCA-GFP/basecalled_data_guppy_301_single_fast5/15/c7f631d9-bdfe-40ab-976b-d2780eb2b63b.fast5"
# Extract read data
read_data <- extract_read_data(file_path = filepath,
                               read_id_fast5_file = NA,
                               plot_debug = FALSE,
                               basecalled_with = 'Guppy',
                               basecall_group = 'Basecall_1D_000',
                               multifast5 = FALSE,
                               model = 'flipflop',
                               plotting_library = 'rbokeh')
fastq <- Biostrings::DNAString(read_data$fastq)

# define parameters for pairwise alignment
match <- 1
mismatch <- -1
type <- 'local'
gapOpening <- 0
gapExtension <- 1
submat <- Biostrings::nucleotideSubstitutionMatrix(match = match,
                                                   mismatch = mismatch)

# define sequences to search for
gfp <- Biostrings::DNAString('GGATCACTCTCGGCATGGACGAGCTGTACAAGTAG')
rc_gfp <- Biostrings::reverseComplement(gfp)

ot <- Biostrings::DNAString("GCCTGTCGCTCTATCTTC")
rc_ot <- Biostrings::DNAString("GAAGATAGAGCGACAGGC")

#ss <- Biostrings::DNAString("AACGATGTACGGCGGG")
#rc_ss <- Biostrings::DNAString("CCCGCCGTACATCGTT")

#rc_ot_ss <- Biostrings::DNAString("GAAGATAGAGCGACAGGCAACGATGTACGGCGGG")
#rc_ss_ot <- Biostrings::DNAString("CCCGCCGTACATCGTTGCCTGTCGCTCTATCTTC")


ot_threshold <- 0.6
rc_ot_threshold <- 0.6
gfp_threshold <- 0.7
rc_gfp_threshold <- 0.7
enable_gfp_search <- FALSE
##########################
### 1. Oligo T search ####
##########################
data_list <- list()
fastq_tmp <- read_data$fastq
# threshold of 0.7 is good for oligoT and RC oligoT
i <- 1
while (TRUE) {
    searchfor <- ot
    as <- Biostrings::pairwiseAlignment(pattern = searchfor,
                                        subject = fastq_tmp,
                                        substitutionMatrix = submat,
                                        type = type,
                                        #scoreOnly = FALSE,
                                        gapOpening = gapOpening,
                                        gapExtension = gapExtension)

    as
    nas <- as@score/searchfor@length

    if (nas < ot_threshold) {
        break
    }

    # substitute NNNs where alignment was found
    aligned_portion <- gsub('-', '', as@subject)
    aligned_portion_length <- nchar(aligned_portion)
    start <- as@subject@range@start
    stop <- start + aligned_portion_length
    substr(fastq_tmp,
           start = start,
           stop = stop) <-
        strrep('N', aligned_portion_length)

    data_list[[i]] <- list(start = start, stop = stop, what_is_it = 'ot')
    i <- i + 1
}
fastq_tmp


#############################
##### 2. RC GFP search ######
#############################
fastq_tmp <- read_data$fastq
while (enable_gfp_search) {
    searchfor <- rc_gfp
    as <- Biostrings::pairwiseAlignment(pattern = searchfor,
                                        subject = fastq_tmp,
                                        substitutionMatrix = submat,
                                        type = type,
                                        #scoreOnly = FALSE,
                                        gapOpening = gapOpening,
                                        gapExtension = gapExtension)

    as
    nas <- as@score/searchfor@length

    if (nas < rc_gfp_threshold) {
        break
    }

    # substitute NNNs where alignment was found
    aligned_portion <- gsub('-', '', as@subject)
    aligned_portion_length <- nchar(aligned_portion)
    start <- as@subject@range@start
    stop <- start + aligned_portion_length
    substr(fastq_tmp,
           start = start,
           stop = stop) <-
        strrep('N', aligned_portion_length)

    data_list[[i]] <- list(start = start, stop = stop, what_is_it = 'rc_gfp')
    i <- i + 1
}
fastq_tmp


#############################
### 3. RC Oligo T search ####
#############################
fastq_tmp <- read_data$fastq
# threshold of 0.7 is good for oligoT and RC oligoT
while (TRUE) {
    searchfor <- rc_ot
    as <- Biostrings::pairwiseAlignment(pattern = searchfor,
                                        subject = fastq_tmp,
                                        substitutionMatrix = submat,
                                        type = type,
                                        #scoreOnly = FALSE,
                                        gapOpening = gapOpening,
                                        gapExtension = gapExtension)

    as
    nas <- as@score/searchfor@length

    if (nas < rc_ot_threshold) {
        break
    }

    # substitute NNNs where alignment was found
    aligned_portion <- gsub('-', '', as@subject)
    aligned_portion_length <- nchar(aligned_portion)
    start <- as@subject@range@start
    stop <- start + aligned_portion_length
    substr(fastq_tmp,
           start = start,
           stop = stop) <-
        strrep('N', aligned_portion_length)

    data_list[[i]] <- list(start = start, stop = stop, what_is_it = 'rc_ot')
    i <- i + 1

}
fastq_tmp


##########################
##### 4. GFP search ######
##########################
fastq_tmp <- read_data$fastq

while (enable_gfp_search) {
    searchfor <- gfp
    as <- Biostrings::pairwiseAlignment(pattern = searchfor,
                                        subject = fastq_tmp,
                                        substitutionMatrix = submat,
                                        type = type,
                                        #scoreOnly = FALSE,
                                        gapOpening = gapOpening,
                                        gapExtension = gapExtension)

    as
    nas <- as@score/searchfor@length

    if (nas < gfp_threshold) {
        break
    }

    # substitute NNNs where alignment was found
    aligned_portion <- gsub('-', '', as@subject)
    aligned_portion_length <- nchar(aligned_portion)
    start <- as@subject@range@start
    stop <- start + aligned_portion_length
    substr(fastq_tmp,
           start = start,
           stop = stop) <-
        strrep('N', aligned_portion_length)

    data_list[[i]] <- list(start = start, stop = stop, what_is_it = 'gfp')
    i <- i + 1
}
fastq_tmp

# TODO: IMPORTANT do this only if data_list is not empty
result <- data_list
result <- purrr::map(result, function(.x) tibble::as_tibble(.x))
result <- dplyr::bind_rows(result, .id = "chunk")
result <- dplyr::select(result, -chunk)
result <- result %>% arrange(start, stop)

### now the algorithm
# good reliable transcript segments are only between
# ot and ot, and rc_ot and rc_ot
# we traverse the data frame what_is_it column and find consecutive repeats
# of ot-ot and rc_ot-rc_ot

what_is_it <- result$what_is_it
starts <- result$start
stops <- result$stop

i <- 1
polyat_seqs <- list()
polyat_seq_count <- 0
while (i < length(what_is_it)) {
    if (what_is_it[i] == 'ot' & what_is_it[i + 1] == 'ot' ) {
        # extract the RC transcript sequence in between the oligo T segments
        start_site <- stops[i]
        stop_site <- starts[i + 1]
        fastq_segment <- fastq[start_site:stop_site]

        # if it is too short a segment, then it is probably some erroneous
        # fusion event, ignor it then
        if (stop_site - start_site > 50) {
            polyat_seq_count <- polyat_seq_count + 1
            polyat_seqs[[polyat_seq_count]] <-
                list(
                    start = start_site,
                    end = stop_site,
                    fastq_segment = as.character(fastq_segment),
                    read_type = 'polyT'
                )
        }
    } else if (what_is_it[i] == 'rc_ot' & what_is_it[i + 1] == 'rc_ot' ) {
        # extract the transcript sequence in between the RC oligo T segments
        start_site <- stops[i]
        stop_site <- starts[i + 1]
        fastq_segment <- fastq[start_site:stop_site]

        # if it is too short a segment, then it is probably some erroneous
        # fusion event, ignor it then
        if (stop_site - start_site > 50) {
            polyat_seq_count <- polyat_seq_count + 1
            polyat_seqs[[polyat_seq_count]] <-
                list(
                    start = start_site,
                    end = stop_site,
                    fastq_segment = as.character(fastq_segment),
                    read_type = 'polyA'
                    )
        }
    }
    i <- i + 1
}

polyat_seqs

# TODO: ensure that you are not doing this on empty
result <- polyat_seqs
result <- purrr::map(result, function(.x) tibble::as_tibble(.x))
result <- dplyr::bind_rows(result, .id = "chunk")
result <- dplyr::select(result, -chunk)
result <- result %>% arrange(start, end)



### Do Clustering
polya_df <- result %>% filter(read_type == 'polyA')
polyt_df <- result %>% filter(read_type == 'polyT')

polya_seqs <- Biostrings::DNAStringSet(polya_df$fastq_segment)
polyt_seqs <- Biostrings::DNAStringSet(polyt_df$fastq_segment)

polya_cluster <- IdClusters(myDistMatrix = NULL,
                            method = "inexact",
                            myXStringSet = polya_seqs,
                            cutoff = 0.9,
                            processors = NULL,
                            verbose = TRUE)
polya_cluster

polya_cluster_list <- list()
for (i in sort(unique(polya_cluster$cluster))) {
    polya_cluster_list[[i]] <- polya_df$fastq_segment[polya_cluster$cluster == i]
}


polyt_cluster <- IdClusters(myDistMatrix = NULL,
                            method = "inexact",
                            myXStringSet = polyt_seqs,
                            cutoff = 0.9,
                            processors = NULL,
                            verbose = TRUE)
polyt_cluster

polyt_cluster_list <- list()
for (i in sort(unique(polyt_cluster$cluster))) {
    polyt_cluster_list[[i]] <- polyt_df$fastq_segment[polyt_cluster$cluster == i]
}
polyt_cluster_list

gapOpening <- 0
gapExtension <- 1


# DO the alignment of all clusters
for (i in seq_along(polya_cluster_list)) {
    alignments <- msa::msa(Biostrings::DNAStringSet(polya_cluster_list[[i]]),
                           method = 'ClustalW',
                           gapOpening = gapOpening,
                           gapExtension = gapExtension)
    consensus <- msa::msaConsensusSequence(alignments, ignoreGaps = FALSE)
    print(gsub('[[:punct:]]', '', consensus))
}

for (i in seq_along(polyt_cluster_list)) {
    alignments <- msa::msa(Biostrings::DNAStringSet(polyt_cluster_list[[i]]),
                           method = 'ClustalW',
                           gapOpening = gapOpening,
                           gapExtension = gapExtension)
    consensus <- msa::msaConsensusSequence(alignments, ignoreGaps = FALSE)
    print(gsub('[[:punct:]]', '', consensus))
}

