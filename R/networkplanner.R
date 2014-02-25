require(sp)
require(igraph)

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
#' @return A NetworkPlan object
#' @export
read_networkplan = function(directory_name) {
    base_dir = normalizePath(directory_name)
    metrics_csv <- read.csv(file.path(base_dir, "metrics-local.csv"), skip=1)
    proj4string <- str_split(readLines(file.path(base_dir, "metrics-local.csv"), n=1), ",")[[1]][[1]]
    
}