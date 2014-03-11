
# assign id, x, y to all points in SpatialLinesDataFrame
get_coord_dataframe <- function(sldf) {
    coord_list <- llply(sldf@lines, function(l) { l@Lines[[1]]@coords })
    coord_matrix <- unique(do.call(rbind, coord_list))
    data.frame(x=coord_matrix[,1], y=coord_matrix[,2]) 
}
    
# get adjacency matrix rep of sldf referencing coord_df ids
# TODO:  Not sure why I get some cycles with this method
get_adjacency_matrix <- function(line_matrix, coord_df) {
    p1 <- as.data.frame(line_matrix[1:nrow(line_matrix),1,1:2])
    p2 <- as.data.frame(line_matrix[1:nrow(line_matrix),2,1:2])
    # need the ids to ensure order after the merge
    p1$pid <- as.numeric(row.names(p1))
    p2$pid <- as.numeric(row.names(p2))
    p1_coord_merge <- merge(p1, coord_df, by.x=c(1,2), by.y=c("x","y"))[,c("pid", "id")]
    p2_coord_merge <- merge(p2, coord_df, by.x=c(1,2), by.y=c("x","y"))[,c("pid", "id")]
    # now get the node ids in the correct order
    p1_ids <- p1_coord_merge[with(p1_coord_merge, order(pid)),"id"] 
    p2_ids <- p2_coord_merge[with(p2_coord_merge, order(pid)),"id"] 
    adj_matrix <- matrix(nrow=nrow(coord_df), ncol=nrow(coord_df))
    adj_matrix[1:nrow(coord_df), 1:nrow(coord_df)] <- 0 
    # must be a better way to set the adjacency matrix incidence cells
    # but for now, use index_matrices (both from/to and to/from)
    index_array1 <- array(c(p1_ids, p2_ids), dim=c(length(p1_ids), 2))
    index_array2 <- array(c(p2_ids, p1_ids), dim=c(length(p1_ids), 2))
    adj_matrix[index_array1] <- 1
    adj_matrix[index_array2] <- 1
    adj_matrix
}
 
