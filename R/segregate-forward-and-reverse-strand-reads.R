segregate_forward_and_reverse_strand_reads <- function(bam_file_path,
                                                       fast5_dir){
    #bam_file_path <- '/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20180516_1429_polya_cdna_shield_run1/alignment_to_genome/aln.bam'
    #fast5_dir <- '/export/valenfs/data/processed_data/MinION/20180516_1429_polya_cdna_shield_run1/basecalled_data'
    # only take the Read ID, strand, mapping quality fields from the BAM file
    what <- c("qname", "strand", "mapq")
    param <- Rsamtools::ScanBamParam(what=what)
    bam <- data.frame(Rsamtools::scanBam(bam_file_path, param=param))
    # remove any rows with NAs in strand or mapping q
    bam <- na.omit(bam)
    # take only those reads that map unabiguously to forward or reverse strands
    bam <- dplyr::filter(bam, strand=='-' | strand=='+')
}





