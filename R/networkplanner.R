library(plyr)
library(ggplot2)
library(rgeos)
library(sp)
library(geosphere)
library(maptools)
library(stringr)

# Prioritized grid & branch indexing
dist_fun <- function(points_1, points_2, 
                     projection_type = proj_var)
{
    if(grepl("units=m", projection_type[1,1])) 
    {
        dist <- sqrt(rowSums((points_1 - points_2)^2))
    }else{
        dist <- distCosine(points_1,points_2,r=6371010) 
    } 
    #     dist <- data.frame(dist , id = data_split[[1]][,"id"])    
    return(dist)
}


prioritized.grid <- function(local_df, shape.file, proj_var = proj4)
{
    dist_fun <- function(points_1, points_2, 
                         projection_type = proj_var)
    {
        if(grepl("units=m", projection_type[1,1])) 
        {
            dist <- sqrt(rowSums((points_1 - points_2)^2))
        }else{
            dist <- distCosine(points_1,points_2,r=6371010) 
        } 
        #     dist <- data.frame(dist , id = data_split[[1]][,"id"])    
        return(dist)
    }
    
    #rename X & Y to long and lat so correspond exactly with proposed grid shape file
    names(local_df)[names(local_df)=="X"] <- "long"
    names(local_df)[names(local_df)=="Y"] <- "lat"
    #subset and truncate local file for most relevant values
    local.grid <- subset(local_df, Metric...System=="grid")
    local.grid <- local.grid[,c("Name",
                                "long",
                                "lat",
                                "Metric...System",
                                "Demographics...Projected.household.count")
                             ]
    
    ##2.0 import proposed grid .shp file -> proposed
    # proposed <- readOGR(folder, "networks-proposed")
    
    
    ## 3.0 Fortify (or flatten) proposed in order to evaluate each vertex or point on shape file lines -> proposed.fortified
    #Fortify shape file to allow for analysis
    proposed.fortified <- fortify(shape.file)
    
    ## 4.0 Merge local.grid along with proposed.fortified vertexe nodes 
    ## ensuring ALL non-matchimg proposed.fortified nodes are kept in the merge -> merged
    merged <- merge(local.grid, proposed.fortified, all.y = T)
    merged$fake.node <- F
    merged[which(is.na(merged$Name)), "fake.node"] <- T
    
    
    ## 8.0 Each "fake.node" will define the start of a unique branch and should be given a "branch.id" from 1 to n
    branch <- (1:sum(as.numeric(merged$fake.node == T)))
    merged$branch <- NA 
    merged$branch[which(merged$fake.node == T)] <- branch
    
    ## 9.0 Find all nodes that match with segments, aka in same "group",  of "fake.nodes"
    ## and label them as "start.nodes" preserving the "branch.id" and "distance" values
    ghost.nodes <- subset(merged, (fake.node == TRUE))
    ghost.nodes <- arrange(ghost.nodes, id)
    
    ### 9.1 determine distance between all nodes
    ##Does not work here as the order between nodes are not consistent so moved to step 11.1 & 12.3 
    
    # ^.
    # ^. Not quite sure what are you trying to do in this step, will revise this part later
    # ^.
    
    ## 10.0 Remove/sub-set all settlements that are "start.nodes" from merged  to new dataframe -> candidate.nodes
    candidate.nodes <- subset(merged, 
                              group %in% ghost.nodes$group 
                              & fake.node==FALSE)
    # use arrange function instead ddply
    candidate.nodes <- arrange(candidate.nodes , id)
    
    #preserve unique branch code originating from existing network 
    candidate.nodes$branch <- ghost.nodes$branch
    #remove candidate nodes from non.candidate nodes dataframe
    
    non.candidate.nodes <- subset(merged, !(id %in% ghost.nodes$id))
    
    ## 11.0 Build new dataframe for ranked settlements -> ranked.settlements
    ## 11.1 determine distance between all nodes
    # data_split <- dlply(non.candidate.nodes, .(order))
    dist <- dist_fun(candidate.nodes[,1:2], ghost.nodes[,1:2], proj_var)
    #     dist <-  distVincentySphere(candidate.nodes[,1:2],
    #                                 ghost.nodes[,1:2],
    #                                 r=6371010) #radius of Earth Network Planner uses
    dist <- data.frame(dist, Name = candidate.nodes[,"Name"])
    
    #reassign distance to candidate.nodes dataframe attributing all length to destination or candidate nodes
    new.candidate.nodes <- merge(candidate.nodes, dist)
    #calculate the unitized MV line required per HH to connect each candidate
    new.candidate.nodes <- mutate(new.candidate.nodes, 
                                  MV.line = dist/Demographics...Projected.household.count)
    
    
    
    
    
    
    branch_identify <- function(candidate_df, non_candidate_df)
    {
        ## 12.0 For 1:nobs merged 
        ## 12.1 select the observation from candidate.nodes that has the minimum grid length/HH connected <- best.candidate
        sequence <- as.integer(1)
        ##settlements <- length(candidate.nodes)+length(unique(non_candidate_df$Name))
        ranked.settlements <- data.frame(NULL)
        #groups <- as.vector(NULL)
        #segments <- as.vector(NULL)
        
        while (dim(candidate_df)[1] > 0) {
            candidate <- subset(candidate_df, MV.line == min(MV.line))[1,]
            candidate <- mutate(candidate, 
                                sequence = sequence)
            sequence <- sequence + 1
            ## 12.2 remove candidate observation from candidate.nodes
            candidate_df <- subset(candidate_df, !(Name %in% candidate$Name))
            ## place candidate in next order of ranked.settlements and assign incremental value for "sequence"
            ranked.settlements <- rbind(ranked.settlements, candidate)
            
            ## 12.3 re-build candidate.nodes by pulling in any settlements that share "group" with candidate observation now in ranked.settlements while assigning it the same "branch.id" as candidate    
            #identify all new lines associated with candidate settlement 
            new.segments <- non_candidate_df$id[which((non_candidate_df$Name == candidate$Name))] #& (non_candidate_df$id != candidate$id))]
            if (length(new.segments)>0) 
            {
                #segments <- c(segments, new.segments)
                #creates 12 variable dataframe such that lines connect to candidate node and settlements are not the candidate node themself
                new.candidate.nodes2 <- as.data.frame(subset(non_candidate_df, (id %in% new.segments) &
                                                                 (Name != candidate$Name)))
                #ensure unique branch labeling scheme is in place showing ancestory 
                v1 <- 1:length(new.segments)
                new.candidate.nodes2$branch <- str_c(candidate$branch,"-",v1)
                #remove non_candidate_df that are now candidates 
                #         non_candidate_df <- subset(non_candidate_df, !(Name %in% candidate$Name) & !(id %in% segments))  
                #         non_candidate_df <- subset(non_candidate_df, !(Name %in% candidate$Name) & !(id %in% new.segments))  
                
                #^. only needed to compare with nodes in new.seg because candidate is in new.segments & new.candidate.node2 is alos in new.seg 
                non_candidate_df <- subset(non_candidate_df, !(id %in% new.segments))  
                #calculate the distance to the new.candidate.nodes                                
                #^. use the new distance function 
                new.dist <- dist_fun(candidate[,c("long", "lat")],
                                     new.candidate.nodes2[,c("long", "lat")],
                                     proj_var)
                
                #                 new.dist <-  distVincentySphere(candidate[,c("long", "lat")],
                #                                                 new.candidate.nodes2[,c("long", "lat")],
                #                                                 r=6371010) #radius of Earth Network Planner uses
                new.dist <- data.frame(dist = new.dist, Name = new.candidate.nodes2[,"Name"])
                #makes compatible 14 variable dataframe for candidate.nodes
                new.candidate.nodes2 <- merge(new.candidate.nodes2, new.dist)    
                ##so here is is an issue where 2 disimilar types are being combined to make a list and not a dataframe
                new.candidate.nodes2 <- mutate(new.candidate.nodes2, 
                                               MV.line = dist/Demographics...Projected.household.count
                )                              
                candidate_df <- rbind(candidate_df, new.candidate.nodes2)
                ## 12.4 Return to step 12.1 until merged dataframe is completely sorted into ranked.settlement 
            }
            
        }
        
        output <- vector("list",2)
        output[[1]] <- ranked.settlements
        output[[2]] <- non_candidate_df
        return(output)
    }
    
    
    ###Starting point of the subnetwork ranking function 
    sub.network.candidate <- function(input_non_candidate_df)
    {
        #####identify the sub-networks
        non_candidate_df <- input_non_candidate_df
        
        ###here is something on top of the whole loop
        sub.network.groups <- vector("list")
        counter <- 1
        
        ##outer loop stars here
        while(length(unique(non_candidate_df$Name)) != 0)
        {
            candidate.groups<- unique(non_candidate_df$Name)
            
            #initialized with one node
            node.names <- candidate.groups[1]
            exit <- F
            
            # Inner loop starts from here
            while(exit == F)
            {
                #pick all the segment id associated with that node
                connected.segments <- unique(non_candidate_df[which(non_candidate_df$Name %in% node.names),"id"])
                if (length(connected.segments) == 0)
                {
                    exit <- T    
                }
                #search all node name associated with the segment
                connected.nodes <- unique(non_candidate_df[which(non_candidate_df$id %in% connected.segments),"Name"])
                #output the names 
                node.names <- unique(c(node.names, connected.nodes))
                #delete all the finded nodes from subnetwork
                non_candidate_df <- subset(non_candidate_df, !(id %in% connected.segments))
                #candidate.groups<- unique(non_candidate_df$Name)
            }
            sub.network.groups[[counter]] <- subset(input_non_candidate_df, Name %in% node.names)
            counter <- counter + 1
        }
        
        #define the function returns the max demo density 
        max.demo <- function(df, demo_col_name="Demographics...Projected.household.count")
        {
            return(df[which(df[,demo_col_name] == max(df[,demo_col_name])),][1,])
        }
        candidates.nodes <- ldply(sub.network.groups, max.demo)
        candidates.nodes <- arrange(candidates.nodes, desc(Demographics...Projected.household.count))
        branch <- paste("S",(1:length(sub.network.groups)), sep="-")
        candidates.nodes$branch <- branch
        
        candidates.nodes <- subset(candidates.nodes, select=c("Name", "id", "branch"))
        
        get.branch <- function(df)
        {
            df$branch <- NULL
            df <- merge(df, candidates.nodes, by = c("Name", "id"), all.x = T)
            new.candidate.nodes <- subset(df, !is.na(branch))
            new.candidate.nodes <- mutate(new.candidate.nodes, 
                                          dist = 0,
                                          MV.line = 0)
            #         non.candidate.nodes <- subset(df, is.na(branch))
            non.candidate.nodes <- df
            output <- branch_identify(new.candidate.nodes, non.candidate.nodes)
            return(output[[1]])
        }
        ranked.sub.networks <- ldply(sub.network.groups, get.branch)
        
        
        
        return(ranked.sub.networks)
    }
    
    output<- branch_identify(new.candidate.nodes, non.candidate.nodes)
    ranked.settlements <- output[[1]]
    ranked.subnetworks <- sub.network.candidate(output[[2]])
    combined.networks <- rbind.fill(ranked.settlements, ranked.subnetworks)
    
    return(combined.networks)
}


# Grid length attribute function
grid.length.attribute <- function(local_df, grid.shape.file, projection_type)
{
    #################$$$$$$$$$$$$$$$$$$$$$
    ### Simple data manipulation
    # change shape data into flat file 
    # Only settlement categorized in "grid" is connected 
    local_df$settlement_id <- rownames(local_df)
    grid_data <- fortify(grid.shape.file)
    settlement <- local_df[,c("X", "Y", "Metric...System", "settlement_id")]
    names(settlement)[1:2] <- c("long", "lat")
    settlement <- subset(settlement, Metric...System == "grid")
    
    # Merge dots with lines
    new_prop_grid <- merge(grid_data,settlement)
    new_prop_grid_2 <- merge(grid_data,settlement, all.x=T)
    single_connections <- new_prop_grid_2[is.na(new_prop_grid_2$settlement_id),"id"]
    
    # Calculate distance for each grid lines using Cosine distance with Great Circle method
    # Assuming the Big Circle is a sphere rather than a ellipses
    data_split <- dlply(grid_data, .(order))
    
    # Introducint the Boolena flag variable proj4 to detect if there is string "units=m" in the 1st line in Metrics-local
    dist <- dist_fun(data_split[[1]][,1:2], data_split[[2]][,1:2], projection_type)
    dist <- data.frame(dist , id = data_split[[1]][,"id"])
    
    #  Double the length of gird lines that has only one connection wiht settlement
    dist[dist$id %in% single_connections,"dist"] <- dist[dist$id %in% single_connections,"dist"] * 2
    
    # Merge grid line length with dots & line data
    new_prop_grid <- merge(new_prop_grid, dist)
    new_prop_grid <- new_prop_grid[order(new_prop_grid$settlement_id), ]
    
    grid_length_attr <- ddply(new_prop_grid, .(settlement_id), summarize, 
                              half_length =  sum(dist)/2)
    
    local_dist <- merge(local_df, grid_length_attr, by="settlement_id", all.x=T)
    local_dist[, "settlement_id"] <- NULL
    return(local_dist)    
}
