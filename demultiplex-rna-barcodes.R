rm(list = ls())
library(hdf5r)
library(Biostrings)
fast5_dir <- '/export/valenfs/data/processed_data/MinION/20180515_1725_polya_direct_rna_shield/polya_estimation/gfp_aligned_reads'
fast5_files_list <- list.files(path=fast5_dir,
                               pattern="\\.fast5$",
                               recursive=TRUE,
                               full.names=TRUE)
df <- data.frame(read_id=character(), filepath=character(), barcode=numeric(), aln_score=numeric())

for (i in 1:length(fast5_files_list)){
    print(i)
    read_path <- fast5_files_list[i]
    f5_obj <- hdf5r::H5File$new(read_path)
    f5_tree <- f5_obj$ls(recursive=TRUE)
    f5_tree <- f5_tree$name
    raw_read_path <- f5_tree[which(f5_tree == 'Raw/Reads') + 2]
    raw_read_path <- f5_tree[which(f5_tree == 'Raw/Reads') + 1]
    ri <- f5_obj[[raw_read_path]]$attr_open('read_id')$read()
    bct <- 'Analyses/Basecall_1D_000/BaseCalled_template'

    fastq <- f5_obj[[bct]]$open('Fastq')$read()

    #get only the fastq sequnce
    fastq <-strsplit(fastq, split = "\n")
    fastq <- fastq[[1]][2]
    fastq_ss <- substr(fastq, start = 1, stop = 100)
    fastq_ss <- gsub("U", "T", fastq_ss)

    ad_10bp <- 'ATCTCTTGCCGTCGCCACC'
    ad_30bp <- 'AGCATTGTAAGTCGCCACC'
    ad_60bp <- 'TAGCCACCAAGTCGCCACC'
    ad_100bp <- 'GAGATTGATGGTCGCCACC'
    ad_150bp <- 'ATAGTCATGGGTCGCCACC'

    submat <- Biostrings::nucleotideSubstitutionMatrix(match = 1,
                                                    mismatch = -1,
                                                       baseOnly = TRUE)

    aln_score_010 <- Biostrings::pairwiseAlignment(pattern=DNAString(ad_10bp),
                                                   subject=DNAString(fastq_ss),
                                                   substitutionMatrix = submat,
                                                   type='local',
                                                   scoreOnly = TRUE)

    aln_score_030 <- Biostrings::pairwiseAlignment(pattern=DNAString(ad_30bp),
                                                   subject=DNAString(fastq_ss),
                                                   substitutionMatrix = submat,
                                                   type='local',
                                                   scoreOnly = TRUE)

    aln_score_060 <- Biostrings::pairwiseAlignment(pattern=DNAString(ad_60bp),
                                                   subject=DNAString(fastq_ss),
                                                   substitutionMatrix = submat,
                                                   type='local',
                                                   scoreOnly = TRUE)

    aln_score_100 <- Biostrings::pairwiseAlignment(pattern=DNAString(ad_100bp),
                                                   subject=DNAString(fastq_ss),
                                                   substitutionMatrix = submat,
                                                   type='local',
                                                   scoreOnly = TRUE)

    aln_score_150 <- Biostrings::pairwiseAlignment(pattern=DNAString(ad_150bp),
                                                   subject=DNAString(fastq_ss),
                                                   substitutionMatrix = submat,
                                                   type='local',
                                                   scoreOnly = TRUE)


     if (aln_score_010 > 13 &&  aln_score_010 > max(aln_score_030, aln_score_060, aln_score_100, aln_score_150)){
         bc <-  10
         as <- aln_score_010
     } else if  (aln_score_030 > 13 &&  aln_score_030 > max(aln_score_010, aln_score_060, aln_score_100, aln_score_150)){
         bc <-  30
         as <- aln_score_030
     } else if (aln_score_060 > 13 &&  aln_score_060 > max(aln_score_030, aln_score_010, aln_score_100, aln_score_150)){
         bc <-  60
         as <- aln_score_060
     } else if (aln_score_100 > 13 &&  aln_score_100 > max(aln_score_030, aln_score_060, aln_score_010, aln_score_150)){
         bc <-  100
         as <- aln_score_100
     } else if (aln_score_150 > 13 &&  aln_score_150 > max(aln_score_030, aln_score_060, aln_score_100, aln_score_010)){
         bc <-  150
         as <- aln_score_150
     } else {
         bc <-  0
         as <- 0

     }

    de <- list(read_id= ri, filepath=read_path, barcode=bc, aln_score=as)
    df <-  rbind(df, de, stringsAsFactors=FALSE)

}
#readRDS(file, refhook = NULL)
data.table::fwrite(df, '/export/valenfs/data/processed_data/MinION/20180515_1725_polya_direct_rna_shield/polya_estimation/filename_barcodes.csv')

