library(tailfinder)
library(ggplot2)
devtools::load_all()

file_path <- system.file(package="tailfinder",
                         "data-raw",
                         "cdna-polya-reads",
                         "0.fast5")

read_data <- extract_read_data(file_path)
norm_data <- z_normalize(read_data$raw_data)

df = data.frame(x=c(1:length(read_data$raw_data)),
                raw_data=read_data$raw_data,
                moves=read_data$moves_sample_wise_vector)

ggplot(data = df, aes(x = x)) +
    geom_line(aes(y = raw_data)) +
    geom_line(aes(y = moves))
