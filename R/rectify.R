#' Convert negative values in the signal to corresponding positive values.
#'
#' @param data a numeric. A numeric vector containing the raw signal.
#'
#' @return A numeric vector containing the modified signal.
#'
rectify <- function(data){
    data <- abs(data)
    return(data)
}
