#' Title
#'
#' @param data
#' @param window_size
#' @param step_size
#'
#' @return
#' @export
#'
#' @examples
rca_make_slope_and_mean_signals <- function(data,
                                            window_size,
                                            step_size) {
    len_data <- length(data)

    spots <- seq(from = 1, to = (len_data - window_size), by = step_size)
    len_spots <- length(spots)
    mean_data <- vector(mode = "numeric", length = len_spots)
    slope <- vector(mode = "numeric", length = len_spots - 1)

    for (i in seq_len(len_spots - 1)) {
        if (i < len_spots - 2) {
            mean_data[i] <- mean(data[spots[i]:(spots[i + 1] - 1)])
            mean_data[i + 1] <- mean(data[spots[i + 1]:(spots[i + 2] - 1)])
            slope[i] = mean_data[i + 1] - mean_data[i]
        } else {
            mean_data[i] <- mean(data[spots[i]:(spots[i + 1] - 1)])
            mean_data[i + 1] <- mean(data[spots[i + 1]:len_data])
            slope[i] = mean_data[i + 1] - mean_data[i]
        }
    }
    return(list(slope = slope,
                mean_data = mean_data))
}
