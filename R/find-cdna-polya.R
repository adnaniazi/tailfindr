#' Find Poly(A) tail in a cDNA read
#'
#' @param file_path Path of the FAST5 file
#'
#' @return A list of Fast5 file data
#' @export
#'
#' @examples
#' find_cdna_polya('path/to/fast5/file')
find_cdna_polya <- function(file_path='a'){
    # file_path <- system.file(package="tailfinder",
    #                          "data-raw",
    #                          "cdna-polya-reads",
    #                          "0.fast5")
    file_path <- '/Users/adnaniazi/Documents/phd/code/tailfinder/data-raw/cdna-polya-reads/6.fast5'

    read_data <- extract_read_data_hdf5r(file_path)
    norm_data <- z_normalize(read_data$raw_data)
    rectified_data <- rectify(norm_data)
    truncated_data <- truncate_spikes(rectified_data, spike_threshold=2)
    truncated_data <- sgolayfilt(truncated_data, p=1, n = 13)

    smoothed_data_1 <- sliding_window('mean', truncated_data, 90, 1)
    smoothed_data_2 <- reverse_sliding_window('mean', truncated_data, 90, 1)
    smoothed_data_3 <- pmin(smoothed_data_1, smoothed_data_2)
    smoothed_data_4 <- rowMeans(cbind(smoothed_data_1, smoothed_data_2))


    df = data.frame(x=c(1:length(read_data$raw_data)),
                    raw_data=read_data$raw_data,
                    norm_data=norm_data,
                    truncated_data=truncated_data,
                    smoothed_data_1=smoothed_data_1,
                    smoothed_data_2=smoothed_data_2,
                    smoothed_data_3=smoothed_data_3,
                    smoothed_data_4=smoothed_data_4,
                    moves=read_data$moves_sample_wise_vector)

    ggplot2::ggplot(data = df, ggplot2::aes(x = x)) +
    ggplot2::geom_line(ggplot2::aes(y = norm_data), color='red') +
    ggplot2::geom_line(ggplot2::aes(y = truncated_data), color='blue') +
    ggplot2::geom_line(ggplot2::aes(y = smoothed_data_1), color='black') +
    ggplot2::geom_line(ggplot2::aes(y = smoothed_data_2), color='green') +
    ggplot2::geom_line(ggplot2::aes(y = smoothed_data_3), color='orange') +
    ggplot2::geom_line(ggplot2::aes(y = smoothed_data_4), color='orange') +
    ggplot2::geom_line(ggplot2::aes(y = moves)) +
    ggplot2::scale_x_continuous(limits = c(length(smoothed_data_1)-ceiling(length(smoothed_data_1)/1),
                                           length(smoothed_data_1)))

    my_data <- xts::xts(df, order.by=Sys.Date()+1:dim(df)[1])
    my_data <- my_data[, 4:9]
    a <- dygraphs::dygraph(my_data, main = "Timeseries plot")
    dygraphs::dyOptions(a, colors = c("pink", "blue", "black", "green", 'purple', 'orange'))
}

