rm(list = ls())
library(hdf5r)
read_path <- '/Users/adnaniazi/Documents/phd/code/tailfinder/data-raw/cdna-polya-reads/1.fast5'
f5_obj <- H5File$new(read_path)
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


