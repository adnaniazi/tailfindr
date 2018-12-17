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
extract_read_data_hdf5r <- function(read_path, plot_debug=FALSE){
    # extract raw data
    f5_obj <- hdf5r::H5File$new(read_path)
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

    # Fastq
    fastq <- f5_obj[[bct]]$open('Fastq')$read()
    # get only the fastq sequnce (ignoring the quality information)
    fastq <-strsplit(fastq, split = "\n")
    fastq <- fastq[[1]][2]

    if (plot_debug) {
    # make a vector of moves interpolated for every sample i.e., make a sample-wise or per-sample vector of moves
        if (event_data$start[1] !=0) {
            moves_sample_wise_vector <- c(rep(NA, event_data$start[1]-1),
                                          rep(event_data$move*0.25+1.5, each=event_data$length[1]),
                                          rep( NA, length(raw_data) - (utils::tail(event_data$start, n=1)+event_data$length[1]-1)))
        } else {
            moves_sample_wise_vector <- c(rep(event_data$move*0.25+1.5, each=event_data$length[1]),
                                          rep( NA, length(raw_data) - (utils::tail(event_data$start, n=1)+event_data$length[1])))
        }
    } else {
        moves_sample_wise_vector <- rep(NA, length(raw_data))
    }

    # create event length data for tail normalization
    l <- event_data$length[1]
    num_events <- length(event_data$start)
    event_length_vector_1 <-  rep(NA, num_events)
    event_length_vector_2 <-  rep(NA, num_events)
    index <- num_events
    length_count <- 1
    divide_by <- 1

    # handle the first row of the event data
    if (event_data$move[1]==1) {
        event_length_vector_1[1] <- 1
    } else {
        event_length_vector_1[index] <- 0.5
        event_length_vector_2[index] <- 0.5
    }
    for(i in (num_events):2) {
        if (event_data$move[i]==2 & event_data$move[i-1]==0) {
            # record the index
            index <- i
            length_count <- 1
            divide_by <- 2
        } else if (event_data$move[i]==1 & event_data$move[i-1]==0) {
            index <- i
            length_count <- 1
            divide_by <- 1
        } else if (event_data$move[i]==0 & (event_data$move[i-1]==1 | event_data$move[i-1]==2)) {
            length_count <- length_count + 1
            # put the previous record
            if (divide_by == 1) {
                event_length_vector_1[index] = length_count
            } else {
                event_length_vector_1[index] = length_count/2
                event_length_vector_2[index] = length_count/2
            }
        } else if ((event_data$move[i]==2 & event_data$move[i-1]==1) | (event_data$move[i]==2 & event_data$move[i-1]==2)) {
            event_length_vector_1[i] <- 0.5
            event_length_vector_2[i] <- 0.5
        } else if ((event_data$move[i]==1 & event_data$move[i-1]==1) | (event_data$move[i]==1 & event_data$move[i-1]==2))  {
            event_length_vector_1[i] <- 1
        } else if  (event_data$move[i]==0 & event_data$move[i-1]==0) {
            length_count <- length_count + 1
        }
    }

    # multiply moves by length of the event (15 for RNA, 5 for DNA)
    event_length_vector_1 <- event_length_vector_1 * l
    event_length_vector_2 <- event_length_vector_2 * l

    event_data <- dplyr::select(event_data, start, length, move, model_state)
    event_data <- cbind(event_data, event_length_vector_1, event_length_vector_2)

    # combine the two vectors
    event_length_vector <- c(event_length_vector_1, event_length_vector_2)

    # reomve NAs
    event_length_vector <- event_length_vector[!is.na(event_length_vector)]

    # median
    samples_per_nt <- median(event_length_vector)

    # drop useless columns in event data
    event_data <- dplyr::select(event_data, start, move, model_state)

    read_data = list(raw_data = raw_data,
                     event_data = event_data,
                     fastq=fastq,
                     moves_sample_wise_vector = moves_sample_wise_vector,
                     fast5_event_length = event_data$length[1],
                     read_id = read_id,
                     read_number = read_number,
                     channel_number = channel_number,
                     start_time = start_time,
                     sampling_rate = sampling_rate,
                     samples_per_nt = samples_per_nt)
    f5_obj$close_all()
    return(read_data)
}
