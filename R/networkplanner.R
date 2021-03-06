require(sp)
require(igraph)
require(rgdal)
require(stringr)
require(plyr)
require(abind)
require(maptools)
require(RCurl)

#' @include util.R
#' @include sequence_models.R 

# Workaround for using igraph in an s4 object slot
setOldClass("igraph")

# The NetworkPlan object
setClass("NetworkPlan", representation(network="igraph", proj="character"))

#' Download scenario from networkplanner into the given directory
#'
#' @param scenario_number scenario_number on http://networkplanner.modilabs.org.
#' @param directory_name path to write the downloaded scenario into. By default,
#'        and if directory_name==NULL, directory_name will be scenario number within
#'        current working directory.
#'        for Windows user you can either use "\\" or "/" as the path separator,
#'        but a single back-slash "\" won't work in R environment as the path construct.
#'        e.g. "C:\\User\\xxx\\" and "User/xxx/" are legit file path.
#' @param username username to login to http://networkplanner.modilabs.org. If NULL,
#'        we assume that the scenario is public.
#' @param password password associated with the previous username. If NULL, 
#'        we assume that the scenario is public. 
#' @param np_url URL of the network planner instance to download scenario from. By default,
#'        it is http://networkplanner.modilabs.org.
#' @export
download_scenario <- function(scenario_number, directory_name=NULL, username=NULL, password=NULL,
                             np_url='http://networkplanner.modilabs.org/') {
    
    ## TODO: Figure out how to handle case that user didnt give login info 
    ## but the SCENARIO happens to be PRIVATE, can't think of a way to validate
    ## until the zip file is downloaded.

    
    # In condition that user didn't give directory_name
    # Use working directory of R session and Scenario number
    # as the folder to save data 
    if (is.null(directory_name)){
        directory_name <- getwd()
        directory_name <- paste(directory_name, scenario_number, sep="/")
    }
    
    # Standardize/ convert to absolute path
    base_dir <- R.utils::getAbsolutePath(normalizePath(directory_name, winslash="/"))
        
    # Create a Boolean flag indicating if the repo if private
    # error handling for only 1 NULL value for user & pass
    private <- (!is.null(username) & !is.null(password))
    if (!sum(!is.null(username), !is.null(password)) %in% c(0,2)){
        stop("You MUST input BOTH username and password if it is a PRIVATE scenario, 
                and leave user and password blank if it is PUBLIC scenario")
    }
    
    # reconstructing url for the zip file
    scenario_addr <- paste("scenarios", 
                           paste(scenario_number, "zip", sep="."), sep="/")
    full_url <- paste(np_url, scenario_addr, sep="")
    
    # Create Curl handle with cookie jar
    my_curl <- getCurlHandle()
    my_curl <- curlSetOpt(cookiejar="",
                          useragent = "Mozilla/5.0",
                          followlocation = TRUE,
                          curl=my_curl)

    # If it is a private repo, then 
    # login network planner and save session into cookie.jar
    if (private == TRUE){
        log_link <- "people/login_"    
        login_url <- paste(np_url, log_link, sep="")
        pars=list(
            username=username,
            password=password)
        postForm(login_url, .params = pars, curl=my_curl)    
    }
    

    # download scenarios to the tmp.zip file in the R session directory
    f <- CFILE("tmp.zip", mode="wb")
    curlPerform(url = full_url, writedata = f@ref, curl=my_curl)
    RCurl::close(f)
    
    # Assume unzip will create the folder, which seem to be a safe assumption 
    # Now unzip the files into base_dir(directory user provided)
    # remove the zipfile once unzipping is finished
    unzip("tmp.zip", exdir = base_dir)
    file.remove("tmp.zip")
}

#' Read scenario from a directory in the filesystem and
#' create it's NetworkPlan.  Converts scenario into an igraph object
#' 
#' @param directory_name absolute or relative path to a directory from which a NetworkPlan is loaded
#' @return A NetworkPlan object
#' @export
read_networkplan <- function(directory_name) {

    base_dir = R.utils::getAbsolutePath(normalizePath(directory_name, winslash="/"))
    
    # read nodes and assign id
    metrics_df <- read.csv(file.path(base_dir, "metrics-local.csv"), skip=1)
    proj4string <- str_extract(readLines(file.path(base_dir, "metrics-local.csv"), n=1), "[+][^,]*")
    
    # check metrics_df for valid X,Y fields
    stopifnot(all(c("X", "Y") %in% names(metrics_df)))

    # read network
    network_shp <- readOGR(dsn=base_dir, layer="networks-proposed")
     
    # make sure projections match (we match metrics and network on coordinates)
    proj_equal <- function(proj4_a, proj4_b) {
       str_extract(proj4_a, "[+]proj[^ ]*")==str_extract(proj4_b, "[+]proj[^ ]*") 
    }
    stopifnot(proj_equal(proj4string, network_shp@proj4string@projargs))
    create_networkplan(metrics_df, network_shp, proj4string)
}


#' Create a NetworkPlan from the components of a scenario
#' 
#' @param metrics_df a dataframe corresponding to the metrics-local output of 
#'        a scenario
#' @param proposed_network_df SpatialLinesDataFrame corresponding to the 
#'        proposed network output of a scenario
#' @param proj4string string representing the projection (in proj4 format)
#' @return A NetworkPlan object
#' @export
create_networkplan <- function(metrics_df, proposed_network_df, proj4string) {

    segments_ids <- decompose_spatial_lines(proposed_network_df)
    network <- create_graph(metrics_df, segments_ids$segments, segments_ids$ids) 
    
    # At this point the network edges should only have the ID attribute from the network shapefile
    # So, assign the rest of the shapefile attributes to the network
    field_names <- names(proposed_network_df@data)
    edge_df <- get.data.frame(network, what="edges")

    # Assumes that ID can be used as an index into the data.frame via ID+1
    edge_df[, field_names] <- proposed_network_df@data[edge_df$ID+1, field_names]
    # is there a better way to merge in all edge fields?  
    for(nm in field_names) {
        vec <- edge_df[,nm]
        # more ugliness (how to deal w/ factors nicely?)
        if(class(vec)=="factor") { vec <- as.character(vec) }
        network[from=edge_df$from, to=edge_df$to, attr=nm] <- vec
    }

    # assign "Network" attributes
    # NOTE:  nid is the id of the relevant metrics-local node
    #        if the vertex doesn't have one it's "fake"
    V(network)$Network..Is.fake <- is.na(V(network)$nid)
    network <- assign_weights(network, weight_field="distance", proj4string)

    # Remove temp vertex attributes
    network <- remove.vertex.attribute(network, "nid") 

    ## Return the NetworkPlan
    new("NetworkPlan", network=network, proj=proj4string)
}

#' Make a NetworkPlan out of a graph and related vertex datafram
#' @param g an igraph
#' @param vertex_df a dataframe of vertices with indices corresponding to the
#'        from/to attributes of the edges of the igraph
#' @param proj the proj4 projection as a string
#' @return A NetworkPlan
#' @export
create_networkplan_from_graph <- function(g, vertex_df, proj) {

    join_vertex_df_igraph <- function(g, vertex_df) {

        vertex_df$id <- 1:nrow(vertex_df)
        col_names <- setdiff(names(vertex_df), "id")
        col_names <- c("id", col_names)
        vertex_df <- vertex_df[,col_names]
        edge_df <- get.data.frame(g, what="edges")
        graph.data.frame(edge_df, directed=F, vertices=vertex_df)
    }     

    network <- join_vertex_df_igraph(g, vertex_df)
    np <- new("NetworkPlan", network=network, proj=proj)
}

#' Take an undirected NetworkPlan and return one that 
#' guarantees that fake nodes belong to disjoint components
#' and that no cycles exist
#' 
#' @param np A NetworkPlan
#' @return A new "cleaned" NetworkPlan object
clean_networkplan <- function(np) {

    # Keep original network and retain the original vertex id
    # so that we can merge back once the directed trees have been created
    original <- np@network
    V(original)$orig_v_id <- 1:length(V(original))

    # ensure that all subgraphs from fake nodes are disjoint and simplify the result
    network <- remove_paths_between_fakes(original)
    network <- simplify(network, remove.loops=TRUE)
    # Take the min span of components in case there are any cycles remaining in them
    # (can occur when scenarios are merged)
    network <- min_span_components(network, np@proj)

    # map the original edge attributes back
    # create the original -> new graph vertex mapping
    orig_new_v_map <- order(V(network)$orig_v_id)
    mapped <- map_edge_attributes(original, network, orig_new_v_map)
    
    # remove temp attrs
    mapped <- remove.vertex.attribute(mapped, "orig_v_id") 
    mapped <- remove.edge.attribute(mapped, "weight") 
    
    np@network <- mapped
    np
}

# Just choose the 1st vertex for now
default_root_selector <- function(graph) {
    as.numeric(V(graph)[1])
}


#' Take an undirected NetworkPlan and make it "directed"
#' such that all components are trees with directed edges
#' to all vertices from the root
#' 
#' @param np A NetworkPlan
#' @return A new "directed" NetworkPlan object
directed_networkplan <- function(np, subnet_root_selector=default_root_selector) {

    # Keep original network and retain the original vertex id
    # so that we can merge back once the directed trees have been created
    original <- np@network
    V(original)$orig_v_id <- 1:length(V(original))

    network <- create_directed_trees(original, root_selector=subnet_root_selector)
 
    # map the original edge attributes back
    # create the original -> new graph vertex mapping
    orig_new_v_map <- order(V(network)$orig_v_id)
    mapped <- map_edge_attributes(original, network, orig_new_v_map)
    
    # remove temp attrs
    mapped <- remove.vertex.attribute(mapped, "orig_v_id") 
    if("vid" %in% list.vertex.attributes(mapped)) { 
        mapped <- remove.vertex.attribute(mapped, "vid") 
    }

    np@network <- mapped
    np
}

#' Default selector for sequence method
#' This simply selects the node with the smallest vertex id
#' 
#' @param node_df node dataframe to apply selection to
#' @return subset of the node_df (MUST BE A DATAFRAME)
default_selector <- function(node_df) {
    subset(node_df, subset=min(node_df$vid)==node_df$vid)
}

#' Default accumulator for accumulate method
#' This calculates the count of all downstream nodes
#' 
#' @param node_df downstream node dataframe
#' @param edge_df downstream edge dataframe
#' @param g the graph
#' @param vid of the current vertex in the graph
#' @return subset of the node_df
default_accumulator <- function(node_df, edge_df, g, vid) {
    data.frame(num_descendents=nrow(node_df))
}

#' Sequence a NetworkPlan via breadth-first search from roots
#' and a selector function (to select from the "frontier" 
#' of SpatialPoints)
#'
#' @param np a NetworkPlan
#' @param selector function that returns which vertex (by id) in the
#'        "frontier" of the search gets selected next based on 
#'        SpatialPointsDataFrame
#' @return A NetworkPlan whose nodes SpatialPointsDataFrame has a sequence 
#'         column and values
#' @export
setGeneric("sequence_plan", function(np, selector=default_selector, validate=F) standardGeneric("sequence_plan"))
setMethod("sequence_plan", signature(np="NetworkPlan"), 
    function(np, selector=default_selector, validate=F) {

        if(validate) {
            stopifnot(can_sequence(np))
        }

        nodes <- get.data.frame(np@network, what="vertices")
        edges <- get.data.frame(np@network, what="edges")
        # start from nodes where Network..Is.root==TRUE since we don't want to
        # sequence "fake" nodes 
        roots <- as.integer(V(np@network)[V(np@network)$Network..Is.root])
        # eliminate the roots that are not part of the network by
        # looking them up in the verts of the edges of the graph
        all_edge_vertices <- c(edges$from, edges$to)
        frontier <- roots[roots %in% all_edge_vertices]
        # get.data.frame returns vertices ordered by vertex id
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
        # only apply to "real" nodes that are ON the network
        real_nodes <- as.integer(V(np@network)[!V(np@network)$Network..Is.fake])
        real_nodes_with_edges <- real_nodes[real_nodes %in% all_edge_vertices]
        num_real_nodes <- length(real_nodes_with_edges)

        stopifnot(length(node_sequence)==num_real_nodes)
        V(np@network)[node_sequence]$Sequence..Far.sighted.sequence <- 1:num_real_nodes
        np
    }
)

#' Accumulate "downstream" nodal values into a summary value
#' for "this" node. 
#'
#' @param np a NetworkPlan
#' @param accumulator function that summarizes or combines all downstream
#'        node values into a new row dataframe attributed to "this" node
#'        the accumulator function takes node and edge dataframes along with 
#'        the graph and vertex id of the current node
#' @return A NetworkPlan whose network vertices and their associated upstream
#'         edges have new attributes/values as defined by accumulator
#' @export
setGeneric("accumulate", function(np, accumulator=default_accumulator, validate=F) standardGeneric("accumulate"))
setMethod("accumulate", signature(np="NetworkPlan"), 
    function(np, accumulator=default_accumulator, validate=F) {

        if(validate) {
            stopifnot(can_sequence(np))
        }

        # get.data.frame returns vertices ordered by vertex id
        nodes_edges <- get.data.frame(np@network, what="both")
        nodes <- nodes_edges$vertices
        edges <- nodes_edges$edges
        nodes$vid <- 1:nrow(nodes)

        # we want to exclude any "fake" nodes (i.e. those nodes that
        # do not have the same attributes as others) if any
        real_nodes <- subset(nodes, !nodes$Network..Is.fake) 
 
        #inner function to "unravel" the dataframe before passing to the accumulator  
        apply_to_down_nodes <- function(df) { 
            
            down_verts <- subcomponent(np@network, df$vid, mode="out")
            down_nodes <- subset(real_nodes, 
                                 vid %in% down_verts)
            down_edges <- subset(edges, from %in% down_verts)
            # now call accumulator callback
            accumulator(down_nodes, down_edges, np@network, df$vid)
        }

        # get new dataframe of accumulator results 
        result_nodes <- ddply(real_nodes, .(vid), apply_to_down_nodes)

        # join back to nodes
        new_nodes <- merge(nodes, result_nodes, by.x="vid", by.y="vid", all.x=TRUE)

        # join to upstream edges
        # Let's NOT do this...user can join nodes to upstream edges via Upstream.ID later on
        # new_edges <- merge(edges, result_nodes, by.x="to", by.y="vid", all.x=TRUE)

        # reorder the edge names (graph.data.frame needs from, to in 1st 2 cols)
        # edge_names <- c(c("from", "to"), setdiff(names(new_edges), c("from", "to")))
        # new_edges <- new_edges[,edge_names]

        # reorder the vertex names (graph.data.frame needs vertex id in 1st col)
        vertex_names <- c(c("vid"), setdiff(names(new_nodes), c("vid")))
        new_vertices <- new_nodes[,vertex_names]

        g <- graph.data.frame(edges, directed=TRUE, new_vertices)
        # get rid of the name attribute as this leads to confusion
        g <- remove.vertex.attribute(g, "name")

        np@network <- g
        # return the networkplan
        np
    }
)

#' Determine whether a NetworkPlan is valid
#'
#' @param np a NetworkPlan
#' @return A boolean indicating whether the NetworkPlan is valid 
#' @export
is_valid_networkplan <- function(np) {
    if(!(class(np)=="NetworkPlan" && class(np@network)=="igraph")) {
        message("Not a proper NetworkPlan object")
        return(FALSE)
    }
    if(length(V(np@network))) {
        v_attrs <- list.vertex.attributes(np@network)
        required_attrs <- c("Network..Is.fake")
        if(sum(v_attrs %in% required_attrs)!=length(required_attrs)) {
            message("NetworkPlan missing required attributes")
            return(FALSE)
        }
    }
    return(TRUE)
}


#' Determine whether a NetworkPlan meets sequencing criteria.
#' This means that there are no paths between "fake nodes" and
#' that there are no cycles.
#'
#' @param np a NetworkPlan
#' @return A boolean indicating whether the NetworkPlan can be sequenced
#' @export
can_sequence = function(np) {

    # Ensure that the graph is directed
    if(!is.directed(np@network)) {
        message("NetworkPlan graph is not directed")
        return(FALSE)
    }

    # Ensure that there are no paths between fake nodes
    fake_vids <- as.numeric(V(np@network)[V(np@network)$Network..Is.fake])
    if (length(fake_vids) > 0) {
        # create all vertex pairs that we need to check for paths
        set_of_pairs <- t(combn(fake_vids, 2))
        # function to get num_paths between pair of vertices (used within apply below)
        get_num_paths <- function(row) {
            edge.connectivity(np@network, row[1], row[2])
        }
        num_paths <- apply(set_of_pairs, 1, get_num_paths)
        if(any(num_paths > 0)) {
            message("NetworkPlan contains paths between fake nodes")
            return(FALSE)
        }
    }
       
    # Now check whether each component is a tree
    # (each node must have one incoming edge except the root node)
    components <- decompose.graph(np@network)
    test_results <- sapply(components, function(g) { vcount(g)==(ecount(g)+1) })
    if(!all(test_results)) {
        message("NetworkPlan components are not all trees")
        return(FALSE)
    }
    return(TRUE)
}


#' Write line shapefile from networkplan object into the given directory
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
write.NetworkPlan <- function(np, directory_name, 
                             nodeFormat='csv', edgeFormat='shp', includeFake=FALSE) {
    
    base_dir <- R.utils::getAbsolutePath(normalizePath(directory_name, winslash="/"))
    
    # subsetting node_df according to includeFake
    output_df <- get.data.frame(np@network, what="vertices")
    if (includeFake == FALSE && 
        ("Network..Is.fake" %in% list.vertex.attributes(np@network))) {
        output_df <- subset(output_df, !Network..Is.fake)
    }
    
    # getting edge SPLDF from NP object
    output_spldf <- get_edge_spldf(np)
    output_spldf@proj4string <- CRS(np@proj)

    # output files based on choice 
    if (nodeFormat == 'csv'){
        csv_dir <- file.path(base_dir, "networkplan-vertices.csv")
        write.csv(output_df, csv_dir, row.names=FALSE)    
    }
    if (edgeFormat =='shp'){
        spldf_dir <- file.path(base_dir, "networkplan-edges")
        writeLinesShape(output_spldf, spldf_dir)    
    }   
    
}

#' Extract spatialpoint and spatialline dataframes for vertices and edges
#'
#' @param np a NetworkPlan
#' @return a list with vertices spatialdataframe and edges spatiallinesdataframe
#' @export
as_spatial_dataframes <- function(np) {

    vertices_spdf <- get.data.frame(np@network, what="vertices")
    coordinates(vertices_spdf) <- ~X + Y
    vertices_spdf@proj4string <- CRS(np@proj)

    edges_spldf <- get_edge_spldf(np)
    edges_spldf@proj4string <- CRS(np@proj)
    
    list(vertices=vertices_spdf, edges=edges_spldf)
}

   
#' Default rollout_functions for sequence_plan_far
default_sequence_model <- list(accumulator=default_accumulator, 
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
#' @param validate whether to check if the networkplan can be sequenced
#' @return A NetworkPlan whose network vertices have a sequence value based
#'         on the rollout_functions.  The edges of the network should also
#'         have the accumulated_field set to the appropriate value.
#' @seealso \code{\link{sequence_plan}}
#' @seealso \code{\link{accumulate}}
#' @export
setGeneric("sequence_plan_far", function(np, sequence_model=default_sequence_model, validate=T) standardGeneric("sequence_plan_far"))
setMethod("sequence_plan_far", signature(np="NetworkPlan"), 
    function(np, sequence_model=default_sequence_model, validate=T) {

        if(!is.directed(np@network)) {
            np <- directed_networkplan(np)
        }

        if(validate) {
            stopifnot(can_sequence(np))
        }

        selector <- sequence_model$selector
        accumulator <- sequence_model$accumulator

        # accumulate values into accumulated_field
        np <- accumulate(np, accumulator=accumulator)
 
        # sequence it via the sequence_selector and return 
        np <- sequence_plan(np, selector=selector)
        
        # now add sequence number to edges
        # sequence number of "to" node should be associated
        # with the edge (i.e. edge gets sequence of downstream node)
        real_nodes <- which(degree(np@network, mode="in")!=0)
        edge_df <- get.data.frame(np@network, what="edges")
        # make sure all edges have a to reference to a "real" node
        stopifnot(nrow(edge_df)==nrow(edge_df[edge_df$to %in% real_nodes,]))
        np@network[from=edge_df$from, to=edge_df$to, 
                   attr="Sequence..Far.sighted.sequence"] <- 
        V(np@network)[edge_df$to]$Sequence..Far.sighted.sequence
        np
    }
)

#' Convert json of form { vertices: {...}, edges: {...} }
#' to an igraph
#' @param json_str json representation of graph as a string
#' @return igraph 
#' @export
json_to_igraph <- function(json_str) {
    js_list <- fromJSON(json_str)
    vdf <- list_of_lists_to_df(js_list$vertices)
    vdf$id <- 1:nrow(vdf)
    col_names <- setdiff(names(vdf), "id")
    col_names <- c("id", col_names)
    vdf <- vdf[,col_names]
    edf <- list_of_lists_to_df(js_list$edges)
    graph.data.frame(edf, FALSE, vdf)
} 
       
#' Convert an igraph to json of form { vertices: {...}, edges: {...} }
#' @param ig igraph
#' @return json representation of an igraph
#' @export
igraph_to_json <- function(ig) {
    vdf <- get.data.frame(ig, what="vertices") 
    edf <- get.data.frame(ig, what="edges") 
    l <- list()
    l$vertices <- df_to_list(vdf)   
    l$edges <- df_to_list(edf)   

    toJSON(l)
}

