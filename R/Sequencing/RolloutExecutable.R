#sort NetworkPlanner's "grid" settlements nodes based on a global priority scheme 
#while denoting each unique branch they're originating from 
#using an implementation of Prim's algorithm on the proposed.grid network.  

setwd("~/github/networkplanner.R/")
source('R/Sequencing/CommonRolloutFunctions.R')

##1.0 - Import metrics.local for only grid-proposed nodes -> local.grid
#load metrics.local to associated settlement points with proposed grid data
local <- read.csv("test_data/metrics-local.csv", skip=1) #RUNTIME ~ 00:28 mins
local$Settlement.id <- rownames(local) #use generic row names for unique ID of each unique settlement point

proj4 <- read.csv("test_data/metrics-local.csv", nrows=1, header = FALSE)

proposed <- readShapeLines("test_data/networks-proposed.shp") #RUNTIME ~ 00:08 mins


#Use output of priortized.grid function as input to far-sighted optimized rollout algorithim 
#takes a shapefile (network) and csv (nodal descriptions and weights) 
#and suggests a sequential, phased roll-out of the system based on a greedy, one step ahead view
#***RUNTIME ~08:00***********
nearsighted_grid <- prioritized.grid.nearsighted(local,proposed, proj4)
##***************************

#Explicitly define greedy grid output as a dataframe
#Sometimes I need to explicitly call the fataframe for greedy.grid - arghhhh
if (length(nearsighted_grid)==2){
  print("Houston, we have a problem with our dataframe")
  nearsighted_grid  <- as.data.frame(nearsighted_grid[1])
}

#Function to determine downstream summations for nearsighted grid AKA Demands and Capacties 
nearsighted_grid_cumulatives <- downstream.sum.calculator(nearsighted_grid)

#Far Sighted function to improve near-sighted greedy grid
#**********Broken!*************
farsighted_grid <- far_sighted_rollout(nearsighted_grid_cumulatives)
#*****

#output csv for analysess 
write.csv(farsighted_grid, "metrics-local-grid-only-rollout_sequence.csv", row.names=F)

# #Output a new shapefile with all the attirbute data of interest
#remove multiple nodes on line segments
metrics_local_with_sequence <- (farsighted_grid[which(!(duplicated(farsighted_grid$id))),])
proposed_with_rollout <- merge(proposed, metrics_local_with_sequence, by.x = "FID", by.y = "id")
writeLinesShape(proposed_with_rollout, "networks-proposed-with-rollout.shp")
