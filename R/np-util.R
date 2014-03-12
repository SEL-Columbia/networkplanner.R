
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
    index_array1 <- array(as.integer(c(p1_ids, p2_ids)), dim=c(length(p1_ids), 2))
    index_array2 <- array(as.integer(c(p2_ids, p1_ids)), dim=c(length(p1_ids), 2))
    adj_matrix[index_array1] <- 1
    adj_matrix[index_array2] <- 1
    adj_matrix
}

get_adjacency_matrix <- function(network_shp){
    line_matrix <- get_coord_matrix(network_shp)
    p1 <- line_matrix[1:dim(line_matrix)[1],1,1:2]
    p2 <- line_matrix[1:dim(line_matrix)[1],2,1:2]
    
    # Merge with node df and detecting ghost nodes & roots
    t1 <- match(data.frame(t(p1)), data.frame(t(nodes@coords)))
    t2 <- match(data.frame(t(p2)), data.frame(t(nodes@coords)))
    match_result <- cbind(t1,t2)
    
    # merge with coordinate matrix(unique) and get the rowname/id
    coord_df<- data.frame(t(unique(rbind(p1,p2))))
    index_array1 <- cbind(match(data.frame(t(p1)), coord_df), 
                     match(data.frame(t(p2)), coord_df))
    index_array2 <- cbind(match(data.frame(t(p2)), coord_df),
                          match(data.frame(t(p1)), coord_df))
    
    
    # Creating the adj matrix
    adj_matrix <- matrix(nrow=ncol(coord_df), ncol=ncol(coord_df), data=0)
    adj_matrix[index_array1] <- 1
    adj_matrix[index_array2] <- 1
    
    return(adj_matrix)
}


my_graph <- graph.adjacency(adj_matrix, mode="undirected")
plot(my_graph)





# Get the co-ordinate matrix out of a SpatialLinesDataFrame (as a three dimensional array)
# Invariant: we should be connecting straight lines in 2D space (ie, 2nd + 3rd dims are 2)
get_coord_matrix = function(sldf) {
    # Looping through every line slot in SLDF object 
    # and save the 2X2 matrix representation of a single edge into 2 3D array. 
    line_coords <- lapply(sldf@lines, function(l) {
        
        lapply(l@Lines, function(segment) {
            my_coords <- segment@coords
            n <- dim(my_coords)[1]
            if(n == 2){
                array(data=my_coords, dim=c(1,2,2))
            }else{
                laply(1:(n-1), function(row_num){
                    cbind(my_coords[row_num, ], my_coords[row_num+1, ])
                })    
            }
            
        })    
    })
    # Looping through every cell in the nested list structure and concatenate the 3D arrays 
    # into the master coordinate array that Adjacency matrix function takes
    line_coords <- do.call(abind, list(lapply(
        line_coords, function(l){ 
            do.call(abind, list(l, along=1))
        }),
        along=1))
    
    # Safety measure to make sure that the coordinate matrix is Nx2x2, 
    # in which N is the number of sigle point in SLDF
    stopifnot(dim(line_coords)[2:3] == c(2,2))
    return(line_coords)
}


# Depricated get_coord_matrix fucntion keeping here for testing purpose for now
# TO DO: turn it into test
# get_coord_matrix = function(sldf) {
#     line_coords <- laply(sldf@lines, function(l) { l@Lines[[1]]@coords })
#     stopifnot(dim(line_coords)[2:3] == c(2,2))
#     line_coords
# }


# To Do: combine the following chunk with adjancey matrix function to increase speed
