reverse_sliding_window <- function(FUN, data, window_size, step_size) {
    #Reverse the data first
    data <- rev(data)

    total <- length(data)
    spots <- seq(from = 1, to = (total - window_size), by = step_size)
    result <- vector(length = length(spots))
    for (i in 1:length(spots)) {
        result[i] <- match.fun(FUN)(data[spots[i]:(spots[i] + window_size - 1)])
    }

    result <- c(result, rep(NA, times=(length(data)-length(result))))
    result <- rev(result)
    return(result)
}
