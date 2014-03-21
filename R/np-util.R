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
    #NOTE:  This vid field is NOT required outside this function
    #       It's only used for accounting between vertices and
    #       segment nodes temporarily
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
    
    #TODO:  do we want directed/undirected here?
    network <- graph.data.frame(edge_df, directed=FALSE, vertex_df)
    #reset names to be consistent with vertex indices
    V(network)$name <- as.numeric(V(network))
    network
}
    
# get 2 subgraphs:
# connected:  The graph reachable from "fake" nodes
# disconnected:  The graph that is not reachable from "fake" nodes
separate_subgraphs <- function(network) {
    fake_vids <- as.numeric(V(network)[is.na(V(network)$nid)])
    reachable <- unique(do.call(c, lapply(fake_vids, function(x) { subcomponent(network, x, mode="ALL") })))
    connected_subgraph <- induced.subgraph(network, reachable)
    unreachable <- setdiff(as.numeric(V(network)), reachable)
    disconnected_subgraph <- induced.subgraph(network, unreachable)
    l <- list(connected=connected_subgraph, 
              disconnected=disconnected_subgraph)
    l
}
 
# Just choose the 1st vertex for now
default_root_selector <- function(graph) {
    as.numeric(V(graph)[1])
}

# From full graph with undirected edges and node ids (i.e. that 
# returned by create_graph), create directed graph with "fake" nodes and
# "root" nodes
create_directed_trees <- function(network, root_selector=default_root_selector) {
    # first, find the connected and disconnected subgraps
    subgraphs <- separate_subgraphs(network)
    connected <- subgraphs$connected
    disconnected <- subgraphs$disconnected

    # now find the roots of the connected subgraph
    fake_vids <- as.numeric(V(connected)[is.na(V(connected)$nid)])

    # Now get all subgraphs from fake vertices for connected net
    # make them directed and union them back together
    get_directed_subgraph <- function(vid, graph) {
        subgraph_vids <- subcomponent(graph, vid, mode="ALL")
        directed <- as.directed(graph, mode="mutual")
        result <- dominator.tree(directed, vid, mode="out")
        dom_tree <- result$domtree
        # I made this quick hack to make it work here
        # TODO: Figure out some solid way to generate directed_subgraph
        directed_subgraph <- induced.subgraph(dom_tree, 
                                              which(!is.na(result$dom)))
        # We should write a test to make sure that 
        # num_vertex(directed_subgraph)  + length(reseult) == num_vertex(directed)
        # and num_vertex(directed_subgraph) == num_non_na_elements(result$dom)
    }
    # Create directed connected subgraph list
    connected_directed_subgraphs <- lapply(fake_vids, 
                                           get_directed_subgraph, 
                                           connected)


    # Temporary visualization for validation within this function
    # TODO: remove in production code base
    # plot(connected_directed_graph, vertex.size=4, edge.arrow.size=0.1)
    
    # Now work on disconnected net
    # Assumes decompose retains vertex ids
    disconnected_subgraphs <- decompose.graph(disconnected)
    disconnected_roots <- sapply(disconnected_subgraphs, root_selector) 
    disconnected_directed_subgraphs <- lapply(disconnected_roots, 
                                              get_directed_subgraph, 
                                              disconnected)
    
    # reduce connected/disconnected directed graphs
    # NOTE:  We do a "disjoint" union here as there will not be
    #        any overlap between graphs (this works much faster with
    #        no messy added attributes)
    #        disjoint runs in seconds rather than minutes for graph.union
    # TODO:  Do we want to keep these separate?  
    connected_directed_graph <- graph.disjoint.union(connected_directed_subgraphs)

    # Designate connected net nodes with is_fake,is_root attributes
    connected_fake_roots <- as.numeric(V(connected_directed_graph)[degree(connected_directed_graph, mode="in")==0])
    connected_real_roots <- as.numeric(V(connected_directed_graph)[nei(connected_fake_roots)])
    V(connected_directed_graph)$is_fake <- ifelse(V(connected_directed_graph) %in% connected_fake_roots,  TRUE, FALSE)
    V(connected_directed_graph)$is_root <- ifelse(V(connected_directed_graph) %in% connected_real_roots,  TRUE, FALSE)

    disconnected_directed_graph <- graph.disjoint.union(disconnected_directed_subgraphs)
    disconnected_real_roots <- as.numeric(V(disconnected_directed_graph)[degree(disconnected_directed_graph, mode="in")==0])
    V(disconnected_directed_graph)$is_fake <- FALSE
    V(disconnected_directed_graph)$is_root <- ifelse(V(disconnected_directed_graph) %in% disconnected_real_roots,  TRUE, FALSE)

    # return the combined directed graph
    combined_directed_graph <- graph.disjoint.union(connected_directed_graph, 
                                                   disconnected_directed_graph)

  
}

# create SPLDF from far sighted sequenced graph with "fake" nodes preserved
get_edge_spldf <- function(np){
    
    # get edge list of vertex representation
    edges <- get.edgelist(np@network)
    edges <- cbind(edge_id=1:nrow(edges), 
                   v1=edges[,1], 
                   v2=edges[,2])
    
    # getting all the edge attributes and assign edge_id based on the sequence
    edge_df <- data.frame(edge.attributes(np@network))
    edge_df$edge_id <- 1:nrow(edge_df)
    
    # creating SpatialLines object from edges matrix 
    splines <- apply(edges, 1, function(row) {
        
        node1 <- as.numeric(row["v1"])
        node2 <- as.numeric(row["v2"])
        edge_id = as.character(row["edge_id"])
        
        Lines(Line(cbind(c(V(x)$X[node1], V(x)$X[node2]),
                         c(V(x)$Y[node1], V(x)$Y[node2]))),
              ID=edge_id)
    })
    
    # creating SpatialLineDataFrame from SpatialLine
    # Because the SpLines are created by the order of each row of edge matrix
    # We can ssume that the order is preserved
    line_df <- SpatialLinesDataFrame(SpatialLines(splines), 
                                     data=data.frame(edge_id=edges[,"edge_id"]))
    line_df <- merge(line_df, edge_df, by="edge_id")
    
    ## A little bit confused about what attributes need to be attached to the SPLDF
    ## Someone can clarify this???
    return(line_df)
}

# Sample benchmarking code
# bench_mark <- microbenchmark(adj1 = get_adjacency_matrix2(network_shp),
#                              adj2 = get_adjacency_matrix(line_matrix=get_coord_matrix(network_shp), coord_df=get_coord_dataframe(network_shp)),
#                              unit="ms", times=1000)
# plot(bench_mark)
# 
# summary(bench_mark)
# qplot(y=time, data=bench_mark, colour=expr) + scale_y_log10()

