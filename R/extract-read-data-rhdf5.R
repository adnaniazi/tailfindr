#' Extract data from ONT Fast5 file
#'
#' Extracts read data from basecalled FAST5 files. The function can't work on
#' raw reads that haven't been basecalled by Albacore.
#'
#' @param read_path Full path of a FAST5 read
#'
#' @return A list of relevant data extracted from the FAST5 file
#'
#' @examples
#' extract_read_data_rhdf5('path/to/fast5/file')

extract_read_data_rhdf5 <- function(read_path){
    # extract raw data
    f5_tree <- rhdf5::h5ls(read_path)
    raw_reads <- f5_tree[(which(f5_tree == "/Raw/Reads") + 1), 1]
    raw_data <- rhdf5::h5read(read_path, raw_reads)$Signal
    ra <- rhdf5::h5readAttributes(read_path, raw_reads)
    ugk <- rhdf5::h5readAttributes(read_path, 'UniqueGlobalKey/channel_id/')

    read_id <- ra$read_id
    read_number <- ra$read_number
    start_time <- ra$start_time

    channel_number <- ugk$channel_number
    sampling_rate <- ugk$sampling_rate

    # extract events data
    tmp_path <- rhdf5::h5ls(read_path)[which(rhdf5::h5ls(read_path) == "/Analyses/Basecall_1D_000/BaseCalled_template")[1], 1]
    event_data <- rhdf5::h5read(read_path, tmp_path)$Events

    # make a vector of moves interpolated for every sample i.e., make a sample-wise or per-sample vector of moves
    if (event_data$start[1] !=0) {
        moves_sample_wise_vector <- c(rep(NA, event_data$start[1]-1),
                                      rep(event_data$move*0.25+1.5, each=event_data$length[1]),
                                      rep( NA, length(raw_data) - (utils::tail(event_data$start, n=1)+event_data$length[1]-1)))
    } else {
        moves_sample_wise_vector <- c(rep(event_data$move*0.25+1.5, each=event_data$length[1]),
                                      rep( NA, length(raw_data) - (utils::tail(event_data$start, n=1)+event_data$length[1])))

    }
    read_data = list(raw_data = raw_data,
                     moves_sample_wise_vector = moves_sample_wise_vector,
                     fast5_event_length = event_data$length[1],
                     read_id = read_id,
                     read_number = read_number,
                     channel_number = channel_number,
                     start_time = start_time,
                     sampling_rate = sampling_rate)
    rhdf5::H5close()
    return(read_data)
}

