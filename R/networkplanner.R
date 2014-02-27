require(sp)
require(igraph)
require(rgdal)

# The NetworkPlan object
setClass("NetworkPlan", representation(nodes="SpatialPointsDataFrame",
                                       network="igraph", existing_network="SpatialLinesDataFrame"))

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
    # TODO: clean up
    coord_matrix <- get_coord_matrix(network_shp)
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