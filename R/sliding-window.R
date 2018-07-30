#' Sliding window
#'
#' Applies a function to window that can slide across a timeseries
#' Acknowledgment: evobiR tool
#'
#' @param FUN Function to apply to each window's data
#' @param data Timeseries data to which to apply the sliding window
#' @param window_size Size of the window
#' @param step_size Step-size of the sliding window
#'
#' @return A signal that has been smoothed using the sliding window. N.B. The
#'   singal is now shorter in lenght
#' @export
#'
#' @examples
#' sliding_window('mean, c(1,2,3,4,5,6,7,8), 2, 1)
#'
sliding_window <- function (FUN, data, window_size, step_size)
{
    total <- length(data)
    spots <- seq(from = 1, to = (total - window_size), by = step_size)
    result <- vector(length = length(spots))
    for (i in 1:length(spots)) {
        result[i] <- match.fun(FUN)(data[spots[i]:(spots[i] + window_size - 1)])
    }
    return(result)
}
