get_fast5_read_id_hdf5r <- function(file_path){
    f5_obj <- hdf5r::H5File$new(file_path)
    f5_tree <- f5_obj$ls(recursive=TRUE)
    f5_tree <- f5_tree$name
    raw_read_path <- f5_tree[which(f5_tree == 'Raw/Reads') + 1]
    print(raw_read_path)
    read_id <- f5_obj[[raw_read_path]]$attr_open('read_id')$read()
    f5_obj$close_all()
    return(list(read_id=read_id, file_path=file_path))
}



