process_optional_params <- function(...){
    if (length(list(...)) > 0) {
        opt_params <- list(...)
        dna_datatype <- opt_params$dna_datatype
    } else {
        dna_datatype <- 'cdna'
    }
    return(dna_datatype)
}

