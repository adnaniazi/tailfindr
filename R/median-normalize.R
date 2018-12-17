#' Title
#'
#' @param timeseries
#'
#' @return
#' @export
#'
#' @examples
median_normalize <- function(timeseries){
    read_median <-  median(timeseries)
    mad <- median(abs(timeseries - read_median))
    median_normalized_data <- (timeseries - read_median)/mad
    return (median_normalized_data)
}



