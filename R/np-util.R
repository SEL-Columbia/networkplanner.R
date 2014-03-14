# assign id, x, y to all points in SpatialLinesDataFrame
get_coord_dataframe <- function(sldf) {
    coord_list <- llply(sldf@lines, function(l) { l@Lines[[1]]@coords })
    coord_matrix <- unique(do.call(rbind, coord_list))
    coord_matrix <- data.frame(x=coord_matrix[,1], y=coord_matrix[,2]) 
    coord_matrix$id <- as.numeric(row.names(coord_matrix))
    coord_matrix
}
    
# get adjacency matrix rep of sldf referencing coord_df ids
# TODO:  Not sure why I get some cycles with this method
get_adjacency_matrix <- function(segment_matrix, segment_node_df) {
    p1 <- as.data.frame(segment_matrix[1:nrow(segment_matrix),1,1:2])
    p2 <- as.data.frame(segment_matrix[1:nrow(segment_matrix),2,1:2])
    # need the ids to ensure order after the merge
    p1$pid <- as.numeric(row.names(p1))
    p2$pid <- as.numeric(row.names(p2))
    p1_coord_merge <- merge(p1, segment_node_df, by.x=c(1,2), by.y=c("X","Y"))[,c("pid", "id")]
    p2_coord_merge <- merge(p2, segment_node_df, by.x=c(1,2), by.y=c("X","Y"))[,c("pid", "id")]
    # now get the node ids in the correct order
    p1_ids <- p1_coord_merge[with(p1_coord_merge, order(pid)),"id"] 
    p2_ids <- p2_coord_merge[with(p2_coord_merge, order(pid)),"id"] 
    adj_matrix <- matrix(nrow=nrow(segment_node_df), ncol=nrow(segment_node_df), data=0)
    
    # must be a better way to set the adjacency matrix incidence cells
    # but for now, use index_matrices (both from/to and to/from)
    
    index_array1 <- cbind(p1_ids, p2_ids)   
    index_array2 <- cbind(p2_ids, p1_ids)

    adj_matrix[index_array1] <- 1
    adj_matrix[index_array2] <- 1
    adj_matrix
}


# Get the co-ordinate matrix out of a SpatialLinesDataFrame (as a three dimensional array)
# Invariant: we should be connecting straight lines in 2D space (ie, 2nd + 3rd dims are 2)
get_segment_matrix = function(sldf) {
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

 
test_edge_pairs <- function(segment_node_df, p1, p2, network) {
    p1_p2 <- data.frame(cbind(p1, p2))
    names(p1_p2) <- c("X1", "Y1", "X2", "Y2")
    p1_p2 <- merge(p1_p2, segment_node_df, by.x=c("X1", "Y1"), by.y=c("X", "Y"), all.x=TRUE)
    p1_p2 <- merge(p1_p2, segment_node_df, by.x=c("X2", "Y2"), by.y=c("X", "Y"), all.x=TRUE)
    edge_df <- get.data.frame(network, what="edges")
    # careful not to test this with large graphs (cartesion product)
    result <- merge(p1_p2, edge_df)
    result <- subset(result, (id.x==from & id.y==to) | (id.x==to & id.y==from))
    expect_equal(nrow(p1_p2)==nrow(result))
    expect_equal(nrow(p1_p2)==(2*nrow(edge_df)))
}
   
# create graph from node_dataframe (NOT sp) and segment_matrix
create_graph <- function(metrics_df, segment_matrix) {
    
    p1 <- segment_matrix[1:dim(segment_matrix)[1],1,1:2]
    p2 <- segment_matrix[1:dim(segment_matrix)[1],2,1:2]
    segment_node_df<- data.frame(unique(rbind(p1,p2)))
    names(segment_node_df) <- c("X", "Y")
    segment_node_df$id <- as.numeric(row.names(segment_node_df))
 
    network_adj_matrix <- get_adjacency_matrix(segment_matrix, segment_node_df)
    network <- graph.adjacency(network_adj_matrix, mode="undirected")
    #vid and vertex index match. add explicit vid field to track
    V(network)$vid <- segment_node_df$id

    # now assign metrics_df attributes back to vertices
    # first, make sure that metrics_df has an nid, representing
    # a unique id for the node (will make it easier to distinguish "fake" nodes)
    metrics_df$nid <- as.numeric(row.names(metrics_df))
    vertex_df <- get.data.frame(network, what=c("vertices"))
    # first get the x,y's to join via
    vertex_df <- merge(vertex_df, segment_node_df, by.x="vid", by.y="id")
    # join to metrics_df and keep all values from vertex_df
    vertex_df <- merge(vertex_df, metrics_df, by.x=c("X","Y"), by.y=c("X","Y"), all.x=TRUE)
    # should be same size (vertices were created from same set of nodes in segment_node_df)
    stopifnot(nrow(vertex_df)==nrow(segment_node_df))
    
    edge_df <- get.data.frame(network, what="edges")
    #reorder vertex df names so that vid is 1st (required to associate with correct edge)
    vertex_names <- c(c("vid"), setdiff(names(vertex_df), c("vid")))
    # vertex_names <- c(c("vid"), names(vertex_df))
    # Because the "vid" column gets lost when it's used as the 1st "linking"
    # column by graph.data.frame, we need to duplicate this column and rename
    # it
    vertex_df <- vertex_df[,vertex_names]
    names(vertex_df)[1] <- "edge_link_id"
    vertex_df$vid <- vertex_df$edge_link_id
    
    #TODO:  do we want directed/undirected here?
    network <- graph.data.frame(edge_df, directed=FALSE, vertex_df)
    network
}
    
#     # Merge with node df and detecting ghost nodes & roots
# #     t1 <- match(data.frame(t(p1)), data.frame(t(nodes@coords)))
# #     t2 <- match(data.frame(t(p2)), data.frame(t(nodes@coords)))
# #     match_result <- cbind(t1,t2)

# Depricated get_coord_matrix fucntion keeping here for testing purpose for now
# TO DO: turn it into test
# get_coord_matrix = function(sldf) {
#     line_coords <- laply(sldf@lines, function(l) { l@Lines[[1]]@coords })
#     stopifnot(dim(line_coords)[2:3] == c(2,2))
#     line_coords
# }


# To Do: combine the following chunk with adjancey matrix function to increase speed
# get_adjacency_matrix2 <- function(network_shp){
#     
#     line_matrix <- get_coord_matrix(network_shp)
#     p1 <- line_matrix[1:dim(line_matrix)[1],1,1:2]
#     p2 <- line_matrix[1:dim(line_matrix)[1],2,1:2]
#     
#     # Merge with node df and detecting ghost nodes & roots
# #     t1 <- match(data.frame(t(p1)), data.frame(t(nodes@coords)))
# #     t2 <- match(data.frame(t(p2)), data.frame(t(nodes@coords)))
# #     match_result <- cbind(t1,t2)
#     
#     # merge with coordinate matrix(unique) and get the rowname/id
# #     coord_df<- data.frame(t(unique(rbind(p1,p2))))
# 
#     coord_df <- t(get_coord_dataframe(network_shp))
#     index_array1 <- cbind(match(data.frame(t(p1)), coord_df), 
#                      match(data.frame(t(p2)), coord_df))
# #     index_array2 <- cbind(match(data.frame(t(p2)), coord_df),
# #                           match(data.frame(t(p1)), coord_df))
#     
#     
#     # Creating the adj matrix
#     adj_matrix <- matrix(nrow=ncol(coord_df), ncol=ncol(coord_df), data=0)
#     adj_matrix[index_array1] <- 1
#     adj_matrix[index_array1[,c(2,1)]] <- 1
#     
#     return(adj_matrix)
# }
# 
# my_graph <- graph.adjacency(get_adjacency_matrix(network_shp), mode="undirected")
# plot(my_graph)
# 
# bench_mark <- microbenchmark(adj1 = get_adjacency_matrix2(network_shp),
#                              adj2 = get_adjacency_matrix(line_matrix=get_coord_matrix(network_shp), coord_df=get_coord_dataframe(network_shp)),
#                              unit="ms", times=1000)
# plot(bench_mark)
# 
# summary(bench_mark)
# qplot(y=time, data=bench_mark, colour=expr) + scale_y_log10()

tmp <- function() { 1 }
