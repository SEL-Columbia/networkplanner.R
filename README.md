networkplanner.R
================

Library to support processing of data associated with [Network Planner](http://networkplanner.modilabs.org).
The core data structure is an [igraph](http://igraph.org)


Installation
---
```
install.packages('devtools')
library(devtools)
install_github("SEL-Columbia/networkplanner.R")
```

Usage Example
---
```
# Set the directory containing the output of a Network Planner
# scenario.  Assumes networkplanner directory is on local machine
base_dir <- "C:/Users/Downloads/2940"

# Generate NetworkPlan object with directed igraph of proposed network
# and nodes from Network Planner
np <- read_networkplan(base_dir)

# Sequence the NetworkPlan object via the 'mv_v_dmd_sequence_model'
np_sequenced <- sequence_plan_far(np, sequence_model=mv_v_dmd_sequence_model)

# Write sequenced NetworkPlan to a directory as nodes (csv) and a network (shp)
write.NetworkPlan(np_sequenced, base_dir)
```

Detailed Overview
---

A NetworkPlan represents a graph-oriented view of the scenario output of 
Network Planner.  The edges of the graph are the segments connecting settlement
nodes (the vertices).  

The call `read_networkplan(base_dir)` returns a NetworkPlan object with
an igraph object in the network slot (i.e. `np@network`).  The igraph returned
is a forest of directed trees.  There are 2 possible roots for these directed
trees:

1.  "Fake" vertices:  

  Vertices that represent the shortest connection from a settlement to an 
  existing network as created by Network Planner.  
  
  These vertices can be found via `V(network)[V(network)$Sequence..Is.fake]`

2.  "Selected" vertices:  

  For trees (subnetworks) that are not connected to an existing network 
  (i.e. have no "Fake" vertex), a root is selected which represents the
  node with maximal demand (we may make this customizable going forward)

Root vertices can be found via `V(network)[V(network)$Sequence..Is.root]`

To create a "Far Sighted Sequence" of the vertices in a NetworkPlan, you can
call `sequence_plan_far(np, sequence_model=mv_v_dmd_sequence_model)`.  This
function takes a custom sequencing algorithm and works in 2 steps:

1.  Accumulate (details in the `accumulate` function in R/networkplanner.R):  

  Given an accumulate function, gather the necessary data associated with each
  vertex (either upstream or downstream as defined by the directionality of the
  graph) and attribute it back to the vertex.  

2.  Sequence (details in the `sequence_plan` function in R/networkplanner.R):  

  Given a selector function, selects the next vertex in the "frontier" of a
  breadth-first-search of the graph in order to generate a sequencing of the
  nodes.   

Data Structures
---

A `NetworkPlan` object consists of the following:
 * `nodes` (a SpatialPointsDataFrame)
 * `network` (an igraph object)
 * `existing_network` (a SpatialLinesDataFrame)

The `nodes` data.frame has the following columns, at least:
 * `id`
 * `is_root`

(note: "fake" nodes from NetworkPlanner are excluded).

Functions
---
 * download_scenario(scenario_number, directory_to_download_into, userpwd=NULL, npURL='http://networkplanner.modilabs.org')
   * scenario_number: scenario number in network planner 
   * directory_to_download_into: a directory to unzip this data into.
   * userpwd: USERNAME:PASSWORD (must be separated by colon). If NULL, scenario_number must be public.
   * npURL: Network Planner URL.
   * _returns_: a `NetworkPlan` object
 * read_networkplan(dirname)
   * dirname: is networkplanner's metrics-local.csv output for a given scenario
   * _returns_: a `NetworkPlan` object
 * sequence_ratio(np, numerator, denominator, nearOrFar='near')
   * np: `NetworkPlan`
   * numerator: fieldName for the numerator (invariant: numerator %in% names(np@nodes))
   * denominator: fieldName for the denominator (invariant: denominator %in% names(np@nodes))
   * nearOrFar: 'near' or 'far'
   * __returns__: `NetworkPlan`, but now there is an additional `sequence` column per node
 * sequence(np, f, nearOrFar='near')
   * np: `NetworkPlan`
   * f: function that takes a data.frame and outputs a numerical objective value
   * nearOrFar: 'near' or 'far'
   * __returns__: `NetworkPlan`, but now there is an additional `sequence` column per node
 * plot.NetworkPlan(np)
   * np: `NetworkPlan`
   * __output__: BOOM! a plot shows up.
 * as.data.frame(np)
   * __returns__: a data frame, where each settlement is represented as a node.
 * as.SpatialLinesDataFrame(np)
   * __returns__: A Spatial Lines Data Frame, where each line has the attributes that are present in the edge attributes of np@network, and all values are NA for lines in the pre-existing network
 * calculateCapacity(np, ??, nearOrFar='near')
   * ??: how do you take values from the node and assign to edges
   * __returns__: A `NetworkPlan`, where np@network now has a `capacity` edge attribute per edge.
 * write.NetworkPlan(np, planDirectory, nodeFormat='csv', edgeFormat='shp')
   * np: `NetworkPlan`
   * planDirectory:  Directory to store nodes and edges within
   
?? branch identify
?? grid-length
?? units

Usage Example
---
```
# assumes networkplanner directory is on local machine
base_dir <- "C:/Users/Downloads/2940"
np <- read_networkplan(base_dir)
# sequence the plan
np_sequenced <- sequence_plan_far(np, sequence_model=mv_v_dmd_sequence_model)
# write sequenced networkplan to a directory 
write.NetworkPlan(np_sequenced, base_dir)
```
