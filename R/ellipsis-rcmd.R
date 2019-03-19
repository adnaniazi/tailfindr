process_optional_params <- function(...){
    # if (length(list(...)) > 0) {
    #     if ("dna_datatype" %in% names(...)) {
    #         dna_datatype <- ...$dna_datatype
    #     } else {
    #         # data parameter is used only when the experiment_type is dna
    #         dna_datatype <- 'cdna'
    #     }
    # } else {
    #     dna_datatype <- 'cdna'
    # }

    if (length(list(...)) > 0) {
        opt_params <- list(...)
        print(opt_params)
        dna_datatype <- opt_params$dna_datatype
    } else {
        dna_datatype <- 'cdna'
    }
    return(dna_datatype)
}

#process_optional_params(a=1, b=1, c=9)
