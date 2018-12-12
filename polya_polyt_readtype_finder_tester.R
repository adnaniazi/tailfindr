library(tailfinder)

polya_folder <- '/export/valenfs/data/processed_data/MinION/20180918_max_cdna_256_polya/gfp_aligned_reads/poly_a_reads'
polyt_folder <- '/export/valenfs/data/processed_data/MinION/20180918_max_cdna_256_polya/gfp_aligned_reads/poly_t_reads/pass/0'

fast5_files_list_polya <- list.files(path=polya_folder,
                                     pattern="\\.fast5$",
                                     recursive=TRUE,
                                     full.names=TRUE)
fast5_files_list_polyt <- list.files(path=polyt_folder,
                                     pattern="\\.fast5$",
                                     recursive=TRUE,
                                     full.names=TRUE)

#dna_tailtype_finder(fast5_files_list_polyt[1])


read_data <- tailfinder::extract_read_data_hdf5r(fast5_files_list_polyt[23], plot_debug=FALSE)
fastq <- read_data$fastq

fa <- Biostrings::DNAString('GGCGTCTGCTTGGGTGTTTAACCTTTTTTTTTTAATGTACTTCGTTCAGTTACGTATTGCT') # GGCGTCTGCTTGGGTGTTTAACCTTTTTTTT
ea <- Biostrings::DNAString('GCAATACGTAACTGAACGAAGT')
# CDNA
fp <- Biostrings::DNAString('TTTCTGTTGGTGCTGATATTGCTGCCATTACGGCCGGG')
ep <- Biostrings::DNAString('GAAGATAGAGCGACAGGCAAGT')

# PCR DNA
#fp <- Biostrings::DNAString('ATTTAGGTGACACTATAGCGCTCCATGCAAACCTGTC')
#ep <- Biostrings::DNAString('CGTTGCCGCCCGGACTC')

rc_fp <- Biostrings::reverseComplement(fp)
rc_ep <- Biostrings::reverseComplement(ep)

match = 1
mismatch = -1
type='local'
gapOpening=0
gapExtension=1

submat <- Biostrings::nucleotideSubstitutionMatrix(match = match,
                                                   mismatch = mismatch,
                                                   baseOnly = TRUE)
adaptors <- list(fa, rc_ep, ea, rc_fp,  fp, ep)


as_fa <- Biostrings::pairwiseAlignment(pattern=fa,
                                    subject=Biostrings::DNAString(substr(fastq, start = 1, stop=100)),
                                    substitutionMatrix = submat,
                                    type=type,
                                    scoreOnly = FALSE,
                                    gapOpening=gapOpening,
                                    gapExtension=gapExtension)

as_fp <- Biostrings::pairwiseAlignment(pattern=fp,
                                       subject=Biostrings::DNAString(substr(fastq, start = 1, stop=100)),
                                       substitutionMatrix = submat,
                                       type=type,
                                       scoreOnly = FALSE,
                                       gapOpening=gapOpening,
                                       gapExtension=gapExtension)


as_ep <- Biostrings::pairwiseAlignment(pattern=ep,
                                       subject=Biostrings::DNAString(substr(fastq, start=nchar(fastq)-50, stop=nchar(fastq))),
                                       substitutionMatrix = submat,
                                       type=type,
                                       scoreOnly = FALSE,
                                       gapOpening=gapOpening,
                                       gapExtension=gapExtension)

as_ea <- Biostrings::pairwiseAlignment(pattern=ea,
                                       subject=Biostrings::DNAString(substr(fastq, start=nchar(fastq)-50, stop=nchar(fastq))),
                                       substitutionMatrix = submat,
                                       type=type,
                                       scoreOnly = FALSE,
                                       gapOpening=gapOpening,
                                       gapExtension=gapExtension)

as_rc_fp <- Biostrings::pairwiseAlignment(pattern=rc_fp,
                                       subject=Biostrings::DNAString(substr(fastq, start=nchar(fastq)-50, stop=nchar(fastq))),
                                       substitutionMatrix = submat,
                                       type=type,
                                       scoreOnly = FALSE,
                                       gapOpening=gapOpening,
                                       gapExtension=gapExtension)

as_rc_ep <- Biostrings::pairwiseAlignment(pattern=rc_ep,
                                       subject=Biostrings::DNAString(substr(fastq, start = 1, stop=70)),
                                       substitutionMatrix = submat,
                                       type=type,
                                       scoreOnly = FALSE,
                                       gapOpening=gapOpening,
                                       gapExtension=gapExtension)

as_fa
as_fp
as_ep
as_ea
as_rc_fp
as_rc_ep

fastq

df <- data.frame(read_data$raw_data)
ggplot2::ggplot(data=df) + ggplot2::geom_line(aes(x=c(1:length(read_data$raw_data)), y=read_data.raw_data))
nchar(fastq)

as_rc_ep@subject@range@start + as_rc_ep@subject@range@width
as_ea@subject@range@start + nchar(fastq)-50
as_ep@subject@range@start + nchar(fastq)-50

#as_rc_ep@subject@range@width
