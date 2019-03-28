#' Explore a fast5 file to find parameters of the experiment
#'
#' This function finds if a read is:
#' \itemize{
#'   \item multifast5 or single fast5
#'   \item 1D read or not-1D
#'   \item basecalled with Albacore or Guppy
#'   \item basecalled using standard model or flipflop model
#'   \item DNA or RNA
#' }
#'
#' @param fast5file_path a character string. Path of a fast5file for determining
#' parameter of the experiment.
#' @param basecall_group a character string. Name of the level
#' in the Fast5 file hierarchy from which to read the data e.g. "Basecall_1D_000"
#'
#' @return A list containing the all the relevant data:
#' \itemize{
#'   \item \code{basecalled_with: 'albacore', 'guppy'}
#'   \item \code{read_is_1d: TRUE, FALSE}
#'   \item \code{model: 'standard', 'flipflop'}
#'   \item \code{fast5type: 'multifast5', 'single'}
#'   \item \code{experiment_type: 'dna', 'rna'}
#' }
#'
#' @examples
#' \dontrun{
#'
#' lst <- explore_basecaller_and_fast5type('/path/to/fast5/file',
#'                                         basecall_group='Basecall_1D_000')
#' }
#'
explore_basecaller_and_fast5type <- function(fast5file_path, basecall_group) {

    f5_obj <- hdf5r::H5File$new(fast5file_path, mode='r')
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
    path_to_check <- paste0('Analyses/', basecall_group)
    if (fast5type == 'single') {
        basecall_1d <- sum(which(f5_tree == path_to_check, arr.ind = T))
    } else {
        basecall_1d <- sum(grepl(path_to_check, f5_tree))
    }
    read_is_1d <- ifelse(basecall_1d > 0, TRUE, FALSE)

    # find if the basecaller is legacy Albacore or the newer Guppy
    if (fast5type == 'multi' & read_is_1d) {
        pth <- paste0('/', path_to_check)
        basecaller_path <- paste(object_name, pth, sep = '')
    } else if (fast5type == 'single' & read_is_1d) {
        basecaller_path <- path_to_check
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
         experiment_type = experiment_type)
}
