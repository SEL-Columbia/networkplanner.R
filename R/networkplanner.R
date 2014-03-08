require(sp)
require(igraph)
require(rgdal)
require(stringr)
require(plyr)
#' @include np-utils.R

# Workaround for using igraph in an s4 object slot
setOldClass("igraph")

# The NetworkPlan object
setClass("NetworkPlan", representation(nodes="SpatialPointsDataFrame", 
                                       network="igraph"))
# Handle case with an existing_network as a subclass for now 
# There's probably a better way to handle it, but without this R complains
# about the existing_network slot being empty if there's nothing there
setClass("NetworkPlanEx", representation(existing_network="SpatialLinesDataFrame"),
                                        contains="NetworkPlan")


#' Download scenario from networkplanner into the given directory
#'
#' @param scenario_number scenario_number on http://networkplanner.modilabs.org
#' @param directory_name path to write the downloaded scenario into. By default,
#'        and if directory_name==NULL, directory_name will be scenario number within
#'        current working directory
#' @param userpwd username and password, separated by : (ex. USER:PASSWORD). If NULL,
#'        we assume that the scenario is public
#' @param np_url URL of the network planner instance to download scenario from. By default,
#'        it is http://networkplanner.modilabs.org
#' @export
download_scenario = function(scenario_number, directory_name=NULL, userpwd=NULL, 
                             np_url='http://networkplanner.modilabs.org/') {
    stop("Not Implemented")
}

#' Read network plan from a directory in the filesystem
#' 
#' @param directory_name absolute or relative path to a directory from which a NetworkPlan is loaded
#' @param debug if TRUE, will verify inputs and run failsafes
#' @return A NetworkPlan object
#' @export
read_networkplan = function(directory_name, debug=F) {
    base_dir = normalizePath(directory_name)
    
    # read nodes and assign id
    metrics_csv <- read.csv(file.path(base_dir, "metrics-local.csv"), skip=1)
    proj4string <- str_extract(readLines(file.path(base_dir, "metrics-local.csv"), n=1), "[+][^,]*")
    nodes <- SpatialPointsDataFrame(coords=metrics_csv[,c("X","Y")], data=metrics_csv,
                                    proj4string=CRS(proj4string))
    nodes$id <- row.names(nodes) # note: this is assumed in the matching code later
    # read network
    network_shp <- readOGR(dsn=base_dir, layer="networks-proposed")
    # TODO: re-project metrics_csv and network_shp to same PROJ? (which one?)
    
    # find roots and separate existing vs. planned network
    
    
    
    p1 <- coord_matrix[1:115,1:2,1]
    p2 <- coord_matrix[1:115,1:2,2]
    t1 <- match(data.frame(t(p1)), data.frame(t(nodes@coords)))
    t2 <- match(data.frame(t(p2)), data.frame(t(nodes@coords)))
    match_result <- cbind(t1,t2)
    
    # TODO: clean up
    line_matrix <- get_coord_matrix(network_shp)
    coord_df <- get_coord_dataframe(network_shp)
    network_adj_matrix <- get_adjacency_matrix(line_matrix, coord_df)
    network_graph <- graph.adjacency(network_adj_matrix, mode="undirected")

    #TODO:  find "roots" and create "dominator trees" from them

    p1 <- as.data.frame(coord_matrix[1:nrow(network_shp),1:2,1])
    p1$FID <- row.names(p1)
    p2 <- as.data.frame(coord_matrix[1:nrow(network_shp),1:2,2])
    p2$FID <- row.names(p2)
    p1 <- merge(p1, as.data.frame(cbind(nodes@coords, id=nodes$id)), 
            by.x=c(1,2), by.y=c("X","Y"), all.x=TRUE)
    p2 <- merge(p2, as.data.frame(cbind(nodes@coords, id=nodes$id)), 
                by.x=c(1,2), by.y=c("X","Y"), all.x=TRUE)
    
    ## TODO: figure out intersection points, with the following objectives:
    ## determine is_root for each node
    ## determine which parts of network_shp go into network::igraph and existing_network::SpatialLinesDataFrame
}

# Get the co-ordinate matrix out of a SpatialLinesDataFrame (as a three dimensional array)
# Invariant: we should be connecting straight lines in 2D space (ie, 2nd + 3rd dims are 2)
get_coord_matrix = function(sldf) {
    line_coords <- laply(sldf@lines, function(l) { t(l@Lines[[1]]@coords) })
    stopifnot(dim(line_coords)[2:3] == c(2,2))
    line_coords
}

#' Default selector for sequence method
#' This simply selects the node with the smallest vertex id
#' 
#' @param node_df SpatialPointsDataFrame to apply selection to
#' @return subset of the node_df
default_selector <- function(node_df) {
    node_df[min(node_df$vid)==node_df$vid,]
}

#' Sequence a NetworkPlan via breadth-first search from roots
#' and a selector function (to select from the "frontier" 
#' of SpatialPoints)
#'
#' @param np a NetworkPlan
#' @param roots the indices of root vertices to sequence from
#' @param selector function that returns which vertex (by id) in the
#'        "frontier" of the search gets selected next based on 
#'        SpatialPointsDataFrame
#' @return A NetworkPlan whose nodes SpatialPointsDataFrame has a sequence 
#'         column and values
#' @export
setGeneric("sequence", function(np, roots, selector=default_selector) standardGeneric("sequence"))
setMethod("sequence", signature(np="NetworkPlan", roots="numeric"), 
    function(np, roots, selector=default_selector) {
        frontier <- roots
        nodes <- np@nodes
        # setup node dataframe with vid (vertex id) field
        # for backrefs
        nodes$vid <- as.numeric(row.names(nodes))
        frontier_df <- nodes[frontier,]
        # keep track of the sequence of the nodes
        selected <- selector(frontier_df)$vid
        # node_sequence a vector whose position represents the sequence index
        # of the node and the value in the position is the node/vertex id
        node_sequence <- selected
        # frontier <- (frontier - selected) + new_neighbors
        frontier <- union(setdiff(frontier, selected), neighbors(np@network, selected))
        while(length(frontier)) {
            frontier_df <- nodes[frontier,]
            # keep track of the sequence of the nodes
            selected <- selector(frontier_df)$vid
            node_sequence <- append(node_sequence, selected)
            # frontier <- (frontier - selected) + new_neighbors
            frontier <- union(setdiff(frontier, selected), neighbors(np@network, selected))
        }
        # now apply the sequence back
        np@nodes[node_sequence, "sequence"] <- 1:length(np@nodes)
        np
    }
)
        
        
        
