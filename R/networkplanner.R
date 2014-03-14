require(sp)
require(igraph)
require(rgdal)
require(stringr)
require(plyr)
require(abind)

#' @include np-utils.R

# Workaround for using igraph in an s4 object slot
setOldClass("igraph")

# The NetworkPlan object
setClass("NetworkPlan", representation(nodes="SpatialPointsDataFrame", 
                                       network="igraph"))
# Handle case with an existing_network as a subclass for now 
# There's probably a better way to handle it, but without this R complains
# about the existing_network slot being empty if there's nothing there
# setClass("NetworkPlanEx", representation(existing_network="SpatialLinesDataFrame"),
#                                        contains="NetworkPlan")


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
    metrics_df <- read.csv(file.path(base_dir, "metrics-local.csv"), skip=1)
    proj4string <- str_extract(readLines(file.path(base_dir, "metrics-local.csv"), n=1), "[+][^,]*")
    nodes <- SpatialPointsDataFrame(coords=metrics_df[,c("X","Y")], data=metrics_df,
                                    proj4string=CRS(proj4string))
    nodes$id <- as.numeric(row.names(nodes)) # note: this is assumed in the matching code later
    
    # read network
    network_shp <- readOGR(dsn=base_dir, layer="networks-proposed")
    # TODO: re-project metrics_csv and network_shp to same PROJ? (which one?)
     
    segment_matrix <- get_segment_matrix(network_shp)
    network <- create_graph(metrics_df, segment_matrix) 

    ## TODO: 
    ## Find all "fake" nodes
    fake_vids <- V(network)[is.na(V(network)$nid)]$vid
    reachable_from_fake <- unique(do.call(c, lapply(fake_vids, function(x) { subcomponent(network, x, mode="ALL") })))
    # reachable_from_fake <- subcomponent(network, fake_vids, mode="ALL")
    unreachable <- setdiff(as.numeric(V(network)), reachable_from_fake)
    unreachable_subgraph <- subgraph(network, unreachable)
      
    ## Create dominator.tree for each fake node
    ## Handle subnets without "fake" nodes
    ## determine is_root for each node
    ## determine which parts of network_shp go into network::igraph and existing_network::SpatialLinesDataFrame
    new("NetworkPlan", nodes=nodes, network=network)
}


#' Default selector for sequence method
#' This simply selects the node with the smallest vertex id
#' 
#' @param node_df SpatialPointsDataFrame to apply selection to
#' @return subset of the node_df
default_selector <- function(node_df) {
    node_df[min(node_df$id)==node_df$id,]
}

#' Default accumulator for accumulate method
#' This calculates the count of all downstream nodes
#' 
#' @param node_df SpatialPointsDataFrame to accumulate values from
#' @return subset of the node_df
default_accumulator <- function(node_df) {
    nrow(node_df)
}

#' Sequence a NetworkPlan via breadth-first search from roots
#' and a selector function (to select from the "frontier" 
#' of SpatialPoints)
#'
#' @param np a NetworkPlan
#'        TODO:  Remove roots once we set them in read method
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
        frontier_df <- nodes[frontier,]
        # keep track of the sequence of the nodes
        selected <- selector(frontier_df)$id
        # node_sequence a vector whose position represents the sequence index
        # of the node and the value in the position is the node/vertex id
        node_sequence <- selected
        # frontier <- (frontier - selected) + new_neighbors
        frontier <- union(setdiff(frontier, selected), neighbors(np@network, selected))
        while(length(frontier)) {
            frontier_df <- nodes[frontier,]
            # keep track of the sequence of the nodes
            selected <- selector(frontier_df)$id
            node_sequence <- append(node_sequence, selected)
            # frontier <- (frontier - selected) + new_neighbors
            frontier <- union(setdiff(frontier, selected), neighbors(np@network, selected))
        }

        # now apply the sequence back
        np@nodes[node_sequence, "sequence"] <- 1:length(np@nodes)
        np
    }
)

#' Accumulate "downstream" nodal values into a summary value
#' for "this" node. 
#' NOTE:  NetworkPlan should be directed from the roots prior
#'        to running. 
#'
#' @param np a NetworkPlan
#'        TODO:  Remove roots once we set them in read method
#' @param roots the indices of root vertices to accumulate from
#' @param accumulated_field name of field to set accumulated value
#' @param accumulator function that summarizes or combines all downstream
#'        node values into a single value attributed to "this" node
#'        (defined by the 'accumulated_field' name)
#'        takes a SpatialPointsDataFrame
#' @return A NetworkPlan whose nodes SpatialPointsDataFrame has an accumulated_field 
#'         column and values
#' @export
setGeneric("accumulate", function(np, roots, accumulated_field, accumulator=default_accumulator) standardGeneric("accumulate"))
setMethod("accumulate", signature(np="NetworkPlan", roots="numeric", accumulated_field="character"), 
    function(np, roots, accumulated_field, accumulator=default_accumulator) {
        #inner function to "unravel" the dataframe before passing to the accumulator  
        apply_to_down_nodes <- function(df) { 
            down_nodes <- data.frame()
            if(length(V(np@network)[df$id])) {
                down_nodes <- np@nodes[subcomponent(np@network, df$id, mode="out"),]
            }
            # now call accumulator callback
            accumulator(down_nodes)
        }
        result <- adply(np@nodes, 1, apply_to_down_nodes)
        # result should be aligned with np@nodes, so just
        # bind the 2nd col of result to the new field
        np@nodes[[accumulated_field]] <- result[,2]
        np
    }
)
