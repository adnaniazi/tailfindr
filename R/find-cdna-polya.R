#' Find Poly(A) tail in a cDNA read
#'
#' @param file_path Path of the FAST5 file
#'
#' @return
#' @export
#'
#' @examples
find_cdna_polya <- function(file_path) {
    # file_path <- system.file(package='tailfinder', 'data-raw', 'cdna-polya-reads', '0.fast5')
    file_path <- "/Users/adnaniazi/Documents/phd/code/tailfinder/data-raw/cdna-polya-reads/0.fast5"
    
    read_data <- extract_read_data(file_path)
    norm_data <- z_normalize(read_data$raw_data)
    
    df = data.frame(x = c(1:length(read_data$raw_data)), raw_data = read_data$raw_data, moves = read_data$moves_sample_wise_vector)
    
    ggplot2::ggplot(data = df, ggplot2::aes(x = x)) + ggplot2::geom_line(ggplot2::aes(y = raw_data)) + ggplot2::geom_line(ggplot2::aes(y = moves))
}
