require(testthat)
require(networkplanner)

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

test_that("reading network plan 174 creates a basically valid NetworkPlan", {
    test_scenario_174_dir <- str_c(test_data_dir, 'Indo_174_Settlements')
    test_scenario_174 <- read_networkplan(test_scenario_dir)
    expect_NetworkPlan_structure(test_scenario_174)
})