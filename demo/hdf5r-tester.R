read_path <- '/Home/ii/adnann/code/tailfinder/data-raw/cdna-polya-reads/107.fast5'

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

# make a vector of moves interpolated for every sample i.e., make a sample-wise or per-sample vector of moves
if (event_data$start[1] !=0) {
    moves_sample_wise_vector <- c(rep(NA, event_data$start[1]-1),
                                  rep(event_data$move*0.25+1.5, each=event_data$length[1]),
                                  rep( NA, length(raw_data) - (utils::tail(event_data$start, n=1)+event_data$length[1]-1)))
} else {
    moves_sample_wise_vector <- c(rep(event_data$move*0.25+1.5, each=event_data$length[1]),
                                  rep( NA, length(raw_data) - (utils::tail(event_data$start, n=1)+event_data$length[1])))

}

# drop useless columns in event data
event_data <- dplyr::select(event_data, start, move, model_state)

read_data = list(raw_data = raw_data,
                 moves_sample_wise_vector = moves_sample_wise_vector,
                 fast5_event_length = event_data$length[1],
                 read_id = read_id,
                 read_number = read_number,
                 channel_number = channel_number,
                 start_time = start_time,
                 sampling_rate = sampling_rate)
f5_obj$close_all()
