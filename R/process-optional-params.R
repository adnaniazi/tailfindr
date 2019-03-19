process_optional_params <- function(...){
    if (length(list(...)) > 0) {
        print('in if')
        opt_params <- list(...)
        print(opt_params)
        dna_datatype <- opt_params$dna_datatype
    } else {
        dna_datatype <- 'cdna'
    }
    return(dna_datatype)
}

