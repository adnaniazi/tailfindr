read_bam <- function(bam_file_path){
    # only take the Read ID, strand, mapping quality fields from the BAM file
    what <- c("qname", "strand", "mapq")
    param <- Rsamtools::ScanBamParam(what=what)
    df_bam <- data.frame(Rsamtools::scanBam(bam_file_path, param=param))
    # remove any rows with NAs in strand or mapping q
    df_bam <- na.omit(df_bam)
    # take only those reads that map unabiguously to forward or reverse strands
    df_bam <- dplyr::filter(df_bam, strand=='-' | strand=='+')
    #c lean up
    df_bam$qname <- as.character(df_bam$qname)
    df_bam$strand <- as.character(df_bam$strand)
    df_bam <- dplyr::rename(df_bam, read_id = qname, bam_mapping_quality=mapq)
    return(df_bam)
}





