rm(list = ls())
library(tailfinder)
files <- list.files('/export/valenfs/data/processed_data/MinION/20181219_max_dna_pcrspikes/renanalysis_with_guppy_flipflop/basecalled_data/workspace/0',
                    recursive = T,
                    full.names = T,
                    pattern = '*.fast5')





extract_read_data_hdf5r_flipflop_multifast5()

rm(list = ls())
library(hdf5r)
library(dplyr)
library(tibble)
read_path_0 <- '/export/valenfs/data/processed_data/MinION/20181219_max_dna_pcrspikes/renanalysis_with_guppy_flipflop/basecalled_data_multifast5/batch_0.fast5'
read_path_1 <- '/export/valenfs/data/processed_data/MinION/20181219_max_dna_pcrspikes/renanalysis_with_guppy_flipflop/basecalled_data_multifast5/batch_1.fast5'
read_path_2 <- '/export/valenfs/data/processed_data/MinION/20181219_max_dna_pcrspikes/renanalysis_with_guppy_flipflop/basecalled_data_multifast5/batch_2.fast5'
read_path_3 <- '/export/valenfs/data/processed_data/MinION/20181219_max_dna_pcrspikes/renanalysis_with_guppy_flipflop/basecalled_data_multifast5/batch_3.fast5'

f5_obj_list <- list()
read_id_obj <- dplyr::tibble(read_id = character(), obj_index = character())

counter <- 1
for (read_path in c(read_path_0, read_path_1, read_path_2, read_path_3)) {
    f5_obj <- hdf5r::H5File$new(read_path, mode='r')
    f5_tree <- f5_obj$ls(recursive=F)
    f5_tree <- f5_tree$name
    f5_tree <- dplyr::mutate(dplyr::tbl_df(f5_tree), obj_index=counter) %>% rename(read_id=value)
    read_id_obj <- rbind(read_id_obj, f5_tree)
    f5_obj_list[[counter]] <- f5_obj
    counter <- counter + 1
}

oi = 2 # object index
ri = 2 # read index
read_id <- read_id_obj[ri,]$read_id
object_idx <- read_id_obj[ri,]$obj_index

# get the object
f5_obj <- f5_obj_list[[object_idx]]

# get the raw data
raw_signal_path <- paste('/', read_id, '/Raw/Signal', sep='')
raw_data <- f5_obj[[raw_signal_path]]$read()

# get the channel number, sampling rate
ch_sr_path <- paste('/', read_id, '/channel_id', sep='')
channel_number <- strtoi(f5_obj[[ch_sr_path]]$attr_open('channel_number')$read())
sampling_rate <- strtoi(f5_obj[[ch_sr_path]]$attr_open('sampling_rate')$read())

# get read_id, read number, start time, and duration
raw_attr_path <-  paste('/', read_id, '/Raw/', sep='')
rid <- f5_obj[[raw_attr_path]]$attr_open('read_id')$read()
read_number <- f5_obj[[raw_attr_path]]$attr_open('read_number')$read()
start_time <- f5_obj[[raw_attr_path]]$attr_open('start_time')$read()
duration <- f5_obj[[raw_attr_path]]$attr_open('duration')$read()

# get eventdata
event_data_fastq_path <- paste('/', read_id, '/Analyses/Basecall_1D_000/BaseCalled_template/', sep='')
event_data <- f5_obj[[event_data_fastq_path]]$open('Events')$read()

# get the start sample of the event table
seg_path <- paste('/', read_id, '/Analyses/Segmentation_000/Summary/segmentation', sep='')
start_sample <- f5_obj[[seg_path]]$attr_open('first_sample_template')$read()

# get the stride and number of events called
basecall_template_path <- paste('/', read_id, '/Analyses/Basecall_1D_000/Summary/basecall_1d_template', sep='')
stride <- f5_obj[[basecall_template_path]]$attr_open('block_stride')$read()
called_events <- f5_obj[[basecall_template_path]]$attr_open('called_events')$read()

# add the start column to the event table for legacy purposes
start_col <-seq(from=start_sample, to=(start_sample + (nrow(event_data)-1)*stride), by=stride)
event_data <- event_data %>% dplyr::mutate(start=start_col)

# get the Fastq sequence
fastq <- f5_obj[[event_data_fastq_path]]$open('Fastq')$read()
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

# median
samples_per_nt <- median(event_length_vector)

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
                 start_sample = start_sample,
                 stride = stride)






