require(testthat)
require(networkplanner)
require(sp)
require(igraph)
require(geosphere)
#' @include np-utils.R

# TO RUN THESE TESTS: setwd(???/networkplanner.R); require(testthat); test_dir('inst/tests')
# TODO: is there a better system?
test_data_dir = "../test_data/"
stopifnot(file.exists(test_data_dir))
proj4str <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84")

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
    df$id <- as.numeric(row.names(df))
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
    roots <- as.numeric(V(np@network)[degree(np@network, mode="in")==0])
    np <- sequence(np, roots)
    e_list <- get.edgelist(np@network)
    from_seq <- np@nodes$sequence[e_list[,1]]
    to_seq <- np@nodes$sequence[e_list[,2]]

    # ensure that all "from" nodes are sequenced before "to" nodes
    expect_equal(length(which(from_seq < to_seq)), length(from_seq))
})

simple_coords_lines <- function() {
    xs <- c(0, 0, 1, 1, 1)
    ys <- c(1, 0, 0, 1, -1)
    xys <- cbind(xs, ys)
    coord_df <- data.frame(x=xs, y=ys)
    coord_df$id <- as.numeric(row.names(coord_df))

    # looks like this:
    # |  |
    #  --
    #    |

    line_matrix <- array(0, dim=c(4, 2, 2))
    line_matrix[1,,] <- xys[1:2,]
    line_matrix[2,,] <- xys[2:3,]
    line_matrix[3,,] <- xys[3:4,]
    line_matrix[4,,] <- xys[c(3,5),]

    list(coord_df=coord_df, line_matrix=line_matrix)
}

simple_NetworkPlan <- function() {
    coords_lines <- simple_coords_lines()
    line_matrix <- coords_lines$line_matrix
    coord_df <- coords_lines$coord_df
    sample_pop <- floor(runif(nrow(coord_df), min=0, max=1000))
    coord_attrs <- data.frame(id=coord_df$id, population=sample_pop)
    sp_df <- SpatialPointsDataFrame(cbind(coord_df$x, coord_df$y), coord_attrs, proj4string=proj4str)
    adj_mat <- get_adjacency_matrix(line_matrix, coord_df)
    network <- graph.adjacency(adj_mat, mode="directed")
    dir_network <- dominator.tree(network, root=1, mode="out")$domtree

    new("NetworkPlan", nodes=sp_df, network=dir_network)
}


test_that("get adjacency matrix is correct", {
     
    coords_lines <- simple_coords_lines()
    line_matrix <- coords_lines$line_matrix
    coord_df <- coords_lines$coord_df
    adj_mat <- get_adjacency_matrix(line_matrix, coord_df)
    
    # create lookup index into adj matrix corresponding to 
    # expected incidence cells
    p1_ids <- c(1, 2, 3, 3)
    p2_ids <- c(2, 3, 4, 5)
    index_array <- array(c(p1_ids, p2_ids), dim=c(length(p1_ids), 2)) 
    expect_equal(adj_mat[index_array], c(1,1,1,1))
})

test_that("accumulator works", {

    np <- simple_NetworkPlan()

    # basic test of default_accumulator
    np <- accumulate(np, 1, "num_descendents")
    expect_equal(np@nodes$num_descendents, c(5, 4, 3, 1, 1))

    # test summing downstream populations 
    sum_pop <- function(df) { sum(df$population) }
    np <- accumulate(np, 1, "sum_pop", accumulator=sum_pop)
    expect_less_than(np@nodes$sum_pop[5], np@nodes$sum_pop[1])

})
