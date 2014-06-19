# read in several scenarios
library(rgeos)

data_dir <- "/home/cjn/np_merge/"
scenarios <- c("bago", "yangon", "ayer")

scen_sp_df_list <- list()
for(scen in scenarios) {
    scen_dir <- paste(data_dir, scen, sep="/")    
    scen_sp_df_list[[scen]] <- readOGR(scen_dir, "networks-proposed")
}
