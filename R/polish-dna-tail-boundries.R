#' Title
#'
#' @param event_data
#' @param read_type
#' @param left_boundary
#' @param right_boundary
#'
#' @return
#' @export
#'
#' @examples
polish_dna_tail_boundries <- function(event_data,
                                      left_boundary,
                                      right_boundary,
                                      read_type){
    if (read_type == 'polyT') {
        row_index_start <- which.min(abs(event_data$start - left_boundary))
        print(row_index_start)

    }

}
