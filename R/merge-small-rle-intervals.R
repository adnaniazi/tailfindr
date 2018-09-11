merge_small_rle_intervals <- function(rle_intersections, threshold){
    #rle_lengths <- c(1123, 8960, 262, 101, 1, 1, 6, 413)
    #rle_values <- c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE)

    rle_lengths <- rle_intersections$lengths
    rle_values <- rle_intersections$values

    len_rle <- length(rle_lengths)

    i <- 0
    new_rle_lengths <- c()
    new_rle_values <- c()
    while (i < len_rle) {
        i <- i + 1
        if (rle_lengths[i] > threshold){
            new_rle_lengths <- c(new_rle_lengths, rle_lengths[i])
            new_rle_values <- c(new_rle_values, rle_values[i])
        }
        else {
            new_rle_lengths[length(new_rle_lengths)] <- new_rle_lengths[length(new_rle_lengths)] + rle_lengths[i]
        }
    }

    # now merge consective same-value intervals
    # $new_rle_lengths
    # [1] 1123 8960  262  109  413
    # $new_rle_values
    # [1]  TRUE FALSE  TRUE FALSE FALSE
    if (length(new_rle_lengths) > 1){
        rle_values <- c(new_rle_values[1])
        rle_lengths <- c(new_rle_lengths[1])
        for (k in 2:length(new_rle_lengths)){
            if (tail(rle_values, n=1) == new_rle_values[k]){
                rle_lengths[rle_lengths==tail(rle_lengths, n=1)] <-
                    new_rle_lengths[k]+tail(rle_lengths, n=1)
            }
            else {
                rle_values <- c(rle_values, new_rle_values[k])
                rle_lengths <- c(rle_lengths, new_rle_lengths[k])
            }
        }
    }

   return( list(lengths=rle_lengths,
                values=rle_values))
}
