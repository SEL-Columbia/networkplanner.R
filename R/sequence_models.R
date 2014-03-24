#' "module" maintaining various downstream accumulator and
#' sequence selector functions for use with the NetworkPlanner.R
#' \code{sequence_plan_far} function
#' 
#' @seealso \code{\link{sequence_plan_far}}
#' @todo incorporate validation check to determine whether model
#'       will work with plan (i.e. does it have the right attributes)

# A simple model computing downstream population and
# sequencing nodes by descending downstream population
#' @export
pop_accumulate <- function(node_df, edge_df, g, vid) { data.frame(sum_pop=sum(df$Pop)) }
#' @export
pop_select_max <- function(df) {
    subset(df, subset=(max(df$sum_pop) == df$sum_pop))
}
#' @export
pop_sequence_model <- list(accumulator=pop_accumulate,
                           selector=pop_select_max)

# Model for computing downstream MV Line per kwh
# and sequencing by selecting the nodes with min values first
#' helper function
mv_v_dmd_get_dist_to_parent <- function(g, vid) {
    distance <- 0
    par <- V(g)[ nei(vid, mode="in") ]
    if(length(par)) {
        distance <- E(g)[ par %->% vid ]$distance
    } 
    distance
}
    
#' @export
mv_v_dmd_accumulate <- function(node_df, edge_df, g, vid) { 
    distance <- mv_v_dmd_get_dist_to_parent(g, vid)
    dmd_yr <- V(g)[vid]$Demand...Projected.nodal.demand.per.year
    all_distances <- c(distance, edge_df$distance)
    sum_distance <- sum(all_distances)
    all_demand <- node_df$Demand...Projected.nodal.demand.per.year
    sum_dmd_yr <- sum(all_demand)
    sum_mv_v_dmd <- -1
    if(all(sum_dmd_yr != 0)) {
        sum_mv_v_dmd <- sum_distance/sum_dmd_yr
    }
    mv_v_dmd <- -1
    if(dmd_yr != 0) {
        mv_v_dmd <- distance/dmd_yr
    } 
    data.frame(distance=distance,
               mv_v_dmd=mv_v_dmd,
               sum_distance=sum_distance, 
               sum_dmd_yr=sum_dmd_yr,
               sum_mv_v_dmd=sum_mv_v_dmd) 
}
#' @export
mv_v_dmd_select_min <- function(df) {
    subset(df, subset=(min(df$sum_mv_v_dmd) == df$sum_mv_v_dmd))
}
#' @export
mv_v_dmd_sequence_model <- list(accumulator=mv_v_dmd_accumulate,
                                selector=mv_v_dmd_select_min)

