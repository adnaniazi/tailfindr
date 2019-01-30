#' Extract data from ONT Fast5 file
#'
#' Extracts read data from basecalled FAST5 files. The function can't work on
#' raw reads that haven't been basecalled by Albacore.
#'
#' @param plot_debug
#' @param read_path Full path of a FAST5 read
#'
#' @return A list of relevant data extracted from the FAST5 file
#' @export
#'
#' @examples
#' extract_read_data_hdf5r('path/to/fast5/file')
extract_read_data_hdf5r_guppy_single_fast5 <- function(read_path,
                                                       plot_debug=FALSE,
                                                       basecalled_with='guppy_standard_model'){
    # extract raw data
    f5_obj <- hdf5r::H5File$new(read_path, mode='r')
    f5_tree <- f5_obj$ls(recursive=TRUE)
    f5_tree <- f5_tree$name
    raw_read_path <- f5_tree[which(f5_tree == 'Raw/Reads') + 2]
    raw_data <- f5_obj[[raw_read_path]]$read()
    ugk <- 'UniqueGlobalKey/channel_id'

    channel_number <- strtoi(f5_obj[[ugk]]$attr_open('channel_number')$read())
    sampling_rate <- strtoi(f5_obj[[ugk]]$attr_open('sampling_rate')$read())

    raw_read_path <- f5_tree[which(f5_tree == 'Raw/Reads') + 1]
    read_id <- f5_obj[[raw_read_path]]$attr_open('read_id')$read()
    read_number <- f5_obj[[raw_read_path]]$attr_open('read_number')$read()
    start_time <- f5_obj[[raw_read_path]]$attr_open('start_time')$read()
    duration <- f5_obj[[raw_read_path]]$attr_open('duration')$read()
    bct <- 'Analyses/Basecall_1D_000/BaseCalled_template'
    event_data <- f5_obj[[bct]]$open('Events')$read()

    seg <- 'Analyses/Segmentation_000/Summary/segmentation'
    start_sample <- f5_obj[[seg]]$attr_open('first_sample_template')$read()

    summary <- 'Analyses/Basecall_1D_000/Summary/basecall_1d_template'
    stride <- f5_obj[[summary]]$attr_open('block_stride')$read()
    called_events <- f5_obj[[summary]]$attr_open('called_events')$read()

    # add the start column to the event table for legacy purposes
    start_col <-seq(from=start_sample, to=(start_sample + (nrow(event_data)-1)*stride), by=stride)
    event_data <- event_data %>% dplyr::mutate(start=start_col)

    # Fastq
    fastq <- f5_obj[[bct]]$open('Fastq')$read()
    # get only the fastq sequnce (ignoring the quality information)
    fastq <-strsplit(fastq, split = "\n")
    fastq <- fastq[[1]][2]

    if (plot_debug) {
    # make a vector of moves interpolated for every sample i.e., make a sample-wise or per-sample vector of moves
        if (start_sample != 0) {
            moves_sample_wise_vector <- c(rep(NA, start_sample-1),
                                          rep(event_data$move*0.25+1.5, each=stride),
                                          rep(NA, length(raw_data) - start_sample - stride*called_events + 1))
        } else {
            moves_sample_wise_vector <- c(rep(event_data$move*0.25+1.5, each=stride),
                                          rep(NA, length(raw_data) - start_sample - stride*called_events + 1))
        }
    } else {
        moves_sample_wise_vector <- rep(NA, length(raw_data))
    }

    # create event length data for tail normalization

    event_length_vector <- rep(NA, called_events)
    count <- 0
    for (i in seq(from=called_events, to=1, by=-1)) {
        if (event_data$move[i] == 1) {
            event_length_vector[i] <- count + 1
            count <- 0
        } else {
            count <- count + 1
        }
    }

    # multiply moves by length of the event (15 for RNA, 5 for DNA)
    event_length_vector <- event_length_vector * stride
    event_data <- cbind(event_data, event_length_vector)

    # reomve NAs
    event_length_vector <- event_length_vector[!is.na(event_length_vector)]

    # Normalizer for flip-flop based data
    samples_per_nt <- mean(event_length_vector[event_length_vector <= quantile(event_length_vector, 0.95)])

    read_data = list(raw_data = raw_data,
                     event_data = event_data,
                     fastq=fastq,
                     moves_sample_wise_vector = moves_sample_wise_vector,
                     fast5_event_length = stride,
                     read_id = read_id,
                     read_number = read_number,
                     channel_number = channel_number,
                     start_time = start_time,
                     sampling_rate = sampling_rate,
                     samples_per_nt = samples_per_nt,
                     start_sample = start_sample)
    f5_obj$close_all()
    return(read_data)
}
