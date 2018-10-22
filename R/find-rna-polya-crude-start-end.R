#' Find the crude start and end of polyA tail in RNA read
#'
#' @return
#'
#' @examples
find_rna_polya_crude_start_end <- function(rle_lengths,
                                           rle_values,
                                           rle_indices,
                                           trough_threshold){

    polya_crude_start <- NA
    polya_crude_end <- NA

    for (i in 1:(length(rle_lengths)-1)) {
        # just check for a big enough trough in the begining of the read followed by a peak in
        # rle instersections
        if (rle_lengths[i] > trough_threshold && !rle_values[i] && rle_values[i+1] )  {
            polya_crude_start <- rle_indices[i]
            polya_crude_end <- rle_indices[i+1]
            break
        }
    }

    return (list(start=polya_crude_start, end=polya_crude_end))
}
