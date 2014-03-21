require(sp)
require(igraph)
require(rgdal)
require(stringr)
require(plyr)
require(abind)
require(maptools)
require(Rcurl)

#' @include np-utils.R

# Workaround for using igraph in an s4 object slot
setOldClass("igraph")

# The NetworkPlan object
setClass("NetworkPlan", representation(network="igraph"))
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
    
    http://networkplanner.modilabs.org/scenarios/1680.html
    scenario_number <- "1680"
    # reconscructing url for the zip file
    scenario_addr <- paste("scenarios", 
                           paste(scenario_number, "zip", sep="."), sep="/")
    full_url <- paste(np_url, scenario_addr, sep="")
    
    # sending request and authenticate
    my_curl <- getCurlHandle(header = TRUE, userpwd = userpwd)
    my_curl <- getCurlHandle()
    tmp_zip_file <- getURL(full_url, curl=my_curl)
    
    
    
    # unzip files 
    f = CFILE("bfile.zip", mode="wb")
    curlPerform(url = full_url, writedata = f@ref)
    close(f)
    
    unzip("bfile.zip", exdir = "./bfile")
    file.remove("bfile.zip")
    
    
    stop("Not yet Implemented")
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
    
    # read network
    network_shp <- readOGR(dsn=base_dir, layer="networks-proposed")
    # TODO: re-project metrics_csv and network_shp to same PROJ? (which one?)
     
    segment_matrix <- get_segment_matrix(network_shp)
    network <- create_graph(metrics_df, segment_matrix) 

    # Now create directed graph from "fake" nodes (for trees connected
    # to the existing network) and "roots" (for trees that are NOT connected)
    # TODO:  Needs testing
    network <- create_directed_trees(network)
    
   
    ## determine which parts of network_shp go into network::igraph and existing_network::SpatialLinesDataFrame
    new("NetworkPlan", network=network)
}


#' Default selector for sequence method
#' This simply selects the node with the smallest vertex id
#' 
#' @param node_df SpatialPointsDataFrame to apply selection to
#' @return subset of the node_df (MUST BE A DATAFRAME)
default_selector <- function(node_df) {
    subset(node_df, subset=min(node_df$vid)==node_df$vid)
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
setGeneric("sequence_plan", function(np, selector=default_selector) standardGeneric("sequence_plan"))
setMethod("sequence_plan", signature(np="NetworkPlan"), 
    function(np, selector=default_selector) {

        # start from nodes where is_root==TRUE since we don't want to
        # sequence "fake" nodes 
        roots <- as.numeric(V(np@network)[V(np@network)$is_root])
        frontier <- roots
        # get.data.frame returns vertices ordered by vertex id
        nodes <- get.data.frame(np@network, what="vertices")
        # setup node dataframe with vid (vertex id) field
        # for backrefs
        nodes$vid <- 1:nrow(nodes)

        # use subset to guarantee a dataframe is returned
        frontier_df <- subset(nodes, subset=nodes$vid %in% frontier)
        # keep track of the sequence of the nodes
        selected <- selector(frontier_df)$vid
        # node_sequence a vector whose position represents the sequence index
        # of the node and the value in the position is the node/vertex id
        node_sequence <- selected
        # frontier <- (frontier - selected) + new_neighbors
        new_neighbors <- as.numeric(V(np@network)[ nei(selected, mode="out") ])
        frontier <- union(setdiff(frontier, selected), new_neighbors)
        while(length(frontier)) {

            frontier_df <- subset(nodes, subset=nodes$vid %in% frontier)
            # keep track of the sequence of the nodes
            selected <- selector(frontier_df)$vid
            node_sequence <- append(node_sequence, selected)
            # frontier <- (frontier - selected) + new_neighbors

            new_neighbors <- as.numeric(V(np@network)[ nei(selected, mode="out") ])
            frontier <- union(setdiff(frontier, selected), new_neighbors)
        }

        # now apply the sequence back
        # only apply to "real" nodes
        num_real_nodes <- length(V(np@network)[!V(np@network)$is_fake])
        stopifnot(length(node_sequence)==num_real_nodes)
        V(np@network)[node_sequence]$sequence <- 1:num_real_nodes
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
#' @param accumulated_field name of field to set accumulated value
#' @param accumulator function that summarizes or combines all downstream
#'        node values into a single value attributed to "this" node
#'        (defined by the 'accumulated_field' name)
#'        takes a SpatialPointsDataFrame
#' @return A NetworkPlan whose nodes SpatialPointsDataFrame has an accumulated_field 
#'         column and values
#' @export
setGeneric("accumulate", function(np, accumulated_field, accumulator=default_accumulator) standardGeneric("accumulate"))
setMethod("accumulate", signature(np="NetworkPlan", accumulated_field="character"), 
    function(np, accumulated_field, accumulator=default_accumulator) {

        # get.data.frame returns vertices ordered by vertex id
        nodes <- get.data.frame(np@network, what="vertices")
        nodes$vid <- 1:nrow(nodes)

        # we want to exclude any "fake" nodes (i.e. those nodes that
        # do not have the same attributes as others) if any
        nodes <- subset(nodes, !nodes$is_fake) 
 
        #inner function to "unravel" the dataframe before passing to the accumulator  
        apply_to_down_nodes <- function(df) { 
            
            down_nodes <- data.frame()
            if(length(V(np@network)[df$vid])) {
                down_nodes <- subset(nodes, 
                                     nodes$vid %in% subcomponent(np@network, df$vid, mode="out"))
            }
            # now call accumulator callback
            accumulator(down_nodes)
        }
        result <- daply(nodes, .(vid), apply_to_down_nodes)

        # Now add the new field back to the vertices
        # result should be aligned with vertices since
        # nodes dataframe is aligned
        g <- set.vertex.attribute(np@network, 
                                  accumulated_field,
                                  index=nodes$vid,
                                  value=result)
        np@network <- g
        # return the networkplan
        np
    }
)


#' Write line shapfile from networkplan object into the given directory
#'
#' @param np a NetworkPlan
#' @param directory_name path to write the downloaded scenario into. By default,
#'        and if directory_name==NULL, directory_name will be scenario number within
#'        current working directory
#' @param nodeFormat a string indicating the type of output file for the nodes in NetworkPlan,
#' only support 'csv' for now.
#' @param edgeFormat a string indicating the type of output file for the edges in NetworkPlan,
#' only support 'shp' for now.
#' @param includeFake a boolean indicating whether to output the fake node in the 
#' vertex/node csv file, Default is set to False
#' @export
write.NetworkPlan = function(np, directory_name, 
                             nodeFormat='csv', edgeFormat='shp', includeFake=FALSE) {
    
    base_dir = normalizePath(directory_name)
    
    # subsetting node_df according to includeFake
    node_df <- get.data.frame(np@network, what="vertices")
    if (includeFake == FALSE){
        output_df <- subset(node_df, !is_fake)
    }
    
    # getting edge SPLDF from NP object
    output_spldf <- get_edge_spldf(np)
    
    # output files based on choice 
    if (nodeFormat == 'csv'){
        csv_dir <- file.path(base_dir, "metrics-local-grid-only-rollout_sequence.csv")
        write.csv(output_df, csv_dir, row.names=FALSE)    
    }
    if (edgeFormat =='shp'){
        spldf_dir <- file.path(base_dir, "metrics-local-grid-only-rollout_sequence.csv")
        writeLinesShape(output_spldf, spldf_dir)    
    }   
    
}
#' Default rollout_functions for sequence_plan_far
default_sequence_model <- list(accumulator=default_accumulator, 
                              accumulated_field="count_downstream",
                              selector=default_selector)


#' Perform a "far-sighted" sequencing by combining calls to
#' accumulate (to accumulate downstream demand for each vertex)
#' and sequence_plan each getting a custom function to perform
#' the accumulation/sequencing
#'
#' @param np a NetworkPlan
#' @param accumulated_field name of field to set accumulated value
#' @param sequence_model a list with 3 members:  
#'        accumulator function to accumulate downstream values (see \code{accumulate})
#'        accumulator_field name of attribute for accumulated value to go in vertices and edges
#'        selector function to perform sequencing (see \code{sequence_plan})
#' @return A NetworkPlan whose network vertices have a sequence value based
#'         on the rollout_functions.  The edges of the network should also
#'         have the accumulated_field set to the appropriate value.
#' @seealso \code{\link{sequence_plan}}
#' @seealso \code{\link{accumulate}}
#' @export
setGeneric("sequence_plan_far", function(np, sequence_model=default_sequence_model) standardGeneric("sequence_plan_far"))
setMethod("sequence_plan_far", signature(np="NetworkPlan", sequence_model="list"), 
    function(np, sequence_model=default_sequence_model) {
        selector <- sequence_model$selector
        accumulator <- sequence_model$accumulator
        accumulated_field <- sequence_model$accumulated_field

        # accumulate values into accumulated_field
        np <- accumulate(np, accumulated_field, accumulator=accumulator)

        # apply the accumulated vertex value to upstream edges
        # only applies to non-root vertices
        non_roots <- as.numeric(V(np@network)[degree(np@network, mode="in")!=0])
        upstream_edges <- as.numeric(E(np@network)[to(non_roots)])
        non_root_values <- get.vertex.attribute(np@network, 
                                                index=non_roots,
                                                accumulated_field)
        g <- set.edge.attribute(np@network, 
                                accumulated_field,
                                index=upstream_edges,
                                value=non_root_values)
        np@network <- g
 
        # sequence it via the sequence_selector and return 
        np <- sequence_plan(np, selector=selector)
        np
    }
)
