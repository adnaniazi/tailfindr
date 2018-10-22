#' Title
#'
#' @param non_polya_region_start
#' @param non_polya_region_end
#' @param moves_sample_wise_vector
#' @param homopolymer_threshold
#'
#' @return
#'
#' @examples
find_non_polya_region_length_without_homopolymers <- function(non_poly_a_seq_start,
                                                              non_poly_a_seq_end,
                                                              moves_sample_wise_vector,
                                                              homopolymer_threshold) {

    data <- moves_sample_wise_vector[non_poly_a_seq_start : non_poly_a_seq_end]
    l1 <- length(data)
    data <- data[!is.na(data)]
    l2 <- length(data)

    compressed.moves_sample_wise_vector <- rep(TRUE, length(data))

    count <- 0
    subtract <- l2-l1

    for (i in 1:length(data)) {
        if (data[i] == 1.5 && count > homopolymer_threshold) {
            compressed.moves_sample_wise_vector[i] <- FALSE
            if (count == homopolymer_threshold+1) {
                subtract <- subtract + homopolymer_threshold
            }
            count <- count + 1
        } else if (data[i] == 1.75 || data[i] == 2) {
            count <- 0
        } else {
            count <- count + 1
        }
    }
    return(length(compressed.moves_sample_wise_vector[compressed.moves_sample_wise_vector==TRUE])-subtract)
}
