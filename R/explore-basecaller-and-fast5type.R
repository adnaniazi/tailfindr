#' Title
#'
#' @param read_path
#'
#' @return
#' @export
#'
#' @examples
explore_basecaller_and_fast5type <- function(read_path) {

    f5_obj <- hdf5r::H5File$new(read_path, mode='r')
    first_read <- f5_obj$open_by_idx(1)
    object_name <- first_read$get_obj_name()

    # find if the fast5 file is single or multifast5 file
    if (grepl('read_', object_name)) {
        fast5type <- 'multi'
        first_read <- f5_obj$open_by_idx(1)
        f5_tree <- first_read$ls(recursive=TRUE)
    } else {
        fast5type <- 'single'
        f5_tree <- f5_obj$ls(recursive=TRUE)
    }
    f5_tree <- f5_tree$name

    # find if the reads are 1D
    if (fast5type == 'single') {
        basecall_1d <- sum(which(f5_tree == 'Analyses/Basecall_1D_000', arr.ind = T))
    } else {
        basecall_1d <- sum(grepl('Analyses/Basecall_1D_000', f5_tree))
    }
    read_is_1d <- ifelse(basecall_1d > 0, TRUE, FALSE)

    # find if the basecaller is legacy Albacore or the newer Guppy
    if (fast5type == 'multi' & read_is_1d) {
        basecaller_path <- paste(object_name, '/Analyses/Basecall_1D_000', sep = '')
    } else if (fast5type == 'single' & read_is_1d) {
        basecaller_path <- 'Analyses/Basecall_1D_000'
    }

    if (read_is_1d) {
        basecaller <- f5_obj[[basecaller_path]]$attr_open('name')$read()
        basecalled_with_albacore <- grepl('Albacore', basecaller)
        basecalled_with_guppy <- grepl('Guppy', basecaller)
        if (basecalled_with_albacore & !basecalled_with_guppy) {
            basecaller <- 'albacore'
        } else if (!basecalled_with_albacore & basecalled_with_guppy) {
            basecaller <- 'guppy'
        }
    } else {
        basecaller <- 'unknown'
    }

    # find if the model used was standard model or flipflop model
    if (basecaller == 'guppy' & read_is_1d) {
        trace <- sum(f5_tree == grepl('.*Trace$', f5_tree), arr.ind = T)
        model <- ifelse(trace > 0, 'flipflop', 'standard')
    } else if (basecaller == 'albacore') {
        model <- 'standard'
    }

    # find if data is dna or rna
    if (fast5type == 'single') {
        context_tags_path <- f5_tree[grepl('.*context_tags$', f5_tree)]
    } else {
        context_tags_path <- paste(object_name,
                            '/context_tags',
                            sep = '')
    }
    sequencing_kit <- f5_obj[[context_tags_path]]$attr_open('sequencing_kit')$read()
    if (grepl('rna', sequencing_kit)) {
        experiment_type <- 'rna'
    } else {
        experiment_type <- 'dna'
    }
    f5_obj$close_all()

    # basecalled with : 'albacore' or 'guppy
    # model : 'standard' or 'flipflop'
    # fast5type: 'multi' or 'single
    # read_is_1d: TRUE or FALSE
    # experiment_type: 'dna' or 'rna'
    list(basecalled_with = basecaller,
         read_is_1d = read_is_1d,
         model = model,
         fast5type = fast5type,
         experiment_type = experiment_type
         )
}
