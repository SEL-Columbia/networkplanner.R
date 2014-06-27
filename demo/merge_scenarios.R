# read in several scenarios
library(rgeos)

data_dir <- "/home/cjn/np_merge/"
scenarios <- c("bago", "yangon", "ayer")

proj_equal <- function(proj4_a, proj4_b) {
   str_extract(proj4_a, "[+]proj[^ ]*")==str_extract(proj4_b, "[+]proj[^ ]*") 
}


scen_sp_df_list <- list()
scen_metrics_df_list <- list()
scen_proj_list <- list()

for(scen in scenarios) {
    scen_dir <- paste(data_dir, scen, sep="/")    
    network_shp <- readOGR(scen_dir, "networks-proposed")
    metrics_df <- read.csv(file.path(scen_dir, "metrics-local.csv"), skip=1)
    proj4string <- str_extract(readLines(file.path(scen_dir, "metrics-local.csv"), n=1), "[+][^,]*")

    # make sure projections match
    stopifnot(proj_equal(proj4string, network_shp@proj4string@projargs))
    scen_proj_list[[scen]] <- proj4string
    scen_sp_df_list[[scen]] <- network_shp
    scen_metrics_df_list[[scen]] <- metrics_df 
}

# make sure projections across scenarios match
stopifnot(length(unique(scen_proj_list))==1)

# merge the scenario network dataframes
union_scen_lines <- Reduce(gUnion, scen_sp_df_list)
union_scen_df <- SpatialLinesDataFrame(union_scen_lines, data.frame(id=length(union_scen_lines@lines[[1]]@Lines)))

# now merge the scenario metrics locals
col_names_per_scen <- lapply(scen_metrics_df_list, names)
shared_col_names <- Reduce(intersect, col_names_per_scen)
outer_merge <- function(df1, df2) { merge(df1, df2, by=shared_col_names, all=T) }
merged_metrics <- Reduce(outer_merge, scen_metrics_df_list)

# order merged set so that we get unique records selecting prioritized duplicates
system_order <- c("grid"=1, "mini-grid"=2, "off-grid"=3, "unelectrified"=4)
merged_metrics <- merged_metrics[order(system_order[merged_metrics$Metric...System]),]
merged_metrics <- merged_metrics[which(!(duplicated(merged_metrics[,c("X","Y")], fromLast=F))),]

np <- create_networkplan(merged_metrics, union_scen_df, scen_proj_list[[1]])

np_clean <- clean_networkplan(np)

np_directed <- directed_networkplan(np, mv_v_dmd_root_selector)

np_sequenced <- sequence_plan_far(np_directed, sequence_model=mv_v_dmd_sequence_model)
