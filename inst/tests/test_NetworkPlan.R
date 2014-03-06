require(testthat)
require(networkplanner)
require(sp)
require(spdep)
require(igraph)

# TO RUN THESE TESTS: setwd(???/networkplanner.R); require(testthat); test_dir('inst/tests')
# TODO: is there a better system?
test_data_dir = "../test_data/"
stopifnot(file.exists(test_data_dir))

# High-level "expect" function, can be used to check the structure of any NetworkPlan
expect_NetworkPlan_structure <- function(np) {
    expect_equal(class(np), "NetworkPlan")
    expect_equal(class(np@network), "igraph")
    expect_equal(class(np@nodes), "SpatialPointsDataFrame")
    expect_equal(class(np@existing_network), "SpatialLinesDataFrame")
    expect_that("id" %in% names(np@nodes))
    expect_that("is_root" %in% names(np@nodes))
}

# temp workaround for creating distance pairs from points
naive_hav_dist_pairs <- function(A, B) {
    result <- matrix(nrow=nrow(A), ncol=nrow(B))
    for (i in 1:nrow(A)) {
        for(j in 1:nrow(B)) {
            result[i, j] <- distHaversine(A[i,], B[j,])
        }
    }
    result
} 

sample_NetworkPlan <- function() { 
    coords <- matrix(runif(20, min=-1.0, max=1.0), nrow=10, ncol=2)
    df <- data.frame(metric=runif(10, min=100000, max=1000000))
    proj4str <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84")
    sp_df <- SpatialPointsDataFrame(coords, df, proj4string=proj4str)
    
    distance_pair_matrix <- naive_hav_dist_pairs(sp_df, sp_df)
    # todo:  construct full graph from distance pairs
    #        then the minimum spanning tree
    #        then add this to the NetworkPlan object
    full_graph <- graph.adjacency(distance_pair_matrix, mode="undirected", weighted=TRUE)
    mst_graph <- minimum.spanning.tree(full_graph)
    
    # create directed tree from mst_graph from an arbitrary root node
    undir_mst <- as.directed(mst_graph, mode="mutual")
    dir_mst_dom <- dominator.tree(undir_mst, root=6, mode="out")
    dir_mst <- dir_mst_dom$domtree

    # There's no existing_network in this plan
    new("NetworkPlan", nodes=sp_df, network=dir_mst)
}

     
    
test_that("reading network plan 174 creates a basically valid NetworkPlan", {
    test_scenario_174_dir <- str_c(test_data_dir, 'Indo_174_Settlements')
    # test_scenario_174 <- read_networkplan(test_scenario_dir)
    # expect_NetworkPlan_structure(test_scenario_174)
})


test_that("sequence of nodes are consistent with graph topology", {
    # net_plan <- read_networkplan(test_scenario_dir)
    # construct NetworkPlan from scratch
    np <- sample_NetworkPlan()
    # strange that a Vertex sequence cannot auto-cast itself to a numeric
    roots <- as.numeric(V(np_s@network)[degree(np_s@network, mode="in")==0])
    np <- sequence(np, roots)
    e_list <- get.edgelist(np@network)
    from_seq <- np@nodes$sequence[e_list[,1]]
    to_seq <- np@nodes$sequence[e_list[,2]]

    # ensure that all "from" nodes are sequenced before "to" nodes
    expect_equal(length(which(from_seq < to_seq)), length(from_seq))
})
    
