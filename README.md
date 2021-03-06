networkplanner.R
================

Library to support processing of data associated with [Network Planner](http://networkplanner.modilabs.org).
The core data structure is an [igraph](http://igraph.org)


Installation
---

NOTE:  This library is no longer compatible with the latest version of igraph.  
The last version of igraph that IS compatible is 0.7.1.

You can install igraph v0.7.1 via the following:
```
install.packages("http://igraph.org/nightly/get/r/igraph_0.7.1.tar.gz", repos=NULL, type="source")
```

Now you can install the rest:

```
install.packages('devtools')
library(devtools)
install_github("SEL-Columbia/networkplanner.R")
```
Windows Users: If you are using a PC, you MUST install Rtools from CRAN, otherwise non-cran packages wont install properly.  Find a free copy of Rtools.exe appropriate for your system below:

<http://cran.r-project.org/bin/windows/Rtools/>

Usage Example
---

Build and sequence NetworkPlan

```
#Load the library first for custom functions
library(networkplanner)

# Set the directory containing the output of a Network Planner
# scenario.  Assumes networkplanner directory is on local machine
base_dir <- "C:/Users/Downloads/2940"

# Generate NetworkPlan object with directed igraph of proposed network
# and nodes from Network Planner
np <- read_networkplan(base_dir)

# Sequence the NetworkPlan object via the 'mv_v_dmd_sequence_model'
np <- sequence_plan_far(np, sequence_model=mv_v_dmd_sequence_model)

# Write sequenced NetworkPlan to a directory as nodes (csv) and a network (shp)
write.NetworkPlan(np, base_dir)
```
Sample plotting code
```
# assumes np is a sequenced NetworkPlan object
# Setup graph for plotting by coloring roots red and labeling
# vertices by their Far.sighted.sequence value
V(np@network)[V(np@network)$Network..Is.root]$color <- "red"
vertex_labels <- get.vertex.attribute(np@network, "Sequence..Far.sighted.sequence")
plot(np@network, vertex.size=4, edge.arrow.size=1, vertex.label=vertex_labels)
```

And the plot
![Sample Plot](http://sel-columbia.github.io/networkplanner.R/img/sample_plot.png)

Notice how the "Sequence roots" in red are not necessarily the roots of their
respective tree/subnetwork.  This is the case when that subnetwork has a "Fake"
root which models the closest connection to an existing grid and is NOT a
settlement node.  

Field Definitions
---

The following are standard fields that are added to a NetworkPlan upon
creation:  

- Vertex Attributes:   
  `Network..Is.fake`:  Whether this vertex is a "Fake" node  

- Edge Attributes:  
  `FID`:  The FID of the corresponding record in the original existing network 
shapefile  
  `distance`:  The distance (in meters) between vertices that this edge spans  

- Directed NetworkPlan Vertex Attributes:
  `Network..Is.root`:  Whether this vertex is a "Sequence root" 

- Sequenced NetworkPlan Vertex Attributes:  
  `Sequence..Far.sighted.sequence`:  The sequence associated with this vertex 
calculated in a "far sighted" manner.  

Custom model fields will vary with the model.  A convention is to prefix
the attribute name with "Sequence".  If this is followed, some sample code
to list these new attributes is:

```
# assumes np is a sequenced NetworkPlan object
vatts <- list.vertex.attributes(np@network) 
vatts[grep("Sequence", vatts)]
```

Detailed Overview
---

A NetworkPlan represents a graph-oriented view of the scenario output of 
Network Planner.  The edges of the graph are the segments connecting settlement
nodes (the vertices).  

The call `read_networkplan(base_dir)` returns a NetworkPlan object with
an igraph object in the network slot (i.e. `np@network`).  The igraph returned
is an undirected graph with vertices representing settlements and edges 
representing connections between.  

In order to sequence this graph, it must be converted into a forest of
directed trees.  You can do this explicitly via `directed_networkplan` or
it will happen implicitly in the call to `sequence_plan_far`.  

Once you have a directed NetworkPlan, there are 2 possible starting vertices 
for each subcomponent of the graph:

1.  "Fake" vertices:  

  Vertices that represent the shortest connection from a settlement to an 
  existing network as created by Network Planner.  
  
  These vertices can be found via `V(network)[V(network)$Network..Is.fake]`

2.  "Selected" vertices:  

  For trees (subcomponents) that are not connected to an existing network 
  (i.e. have no "Fake" vertex), a root is selected which represents the
  node with maximal demand (or some custom root selection function)

"Network root" vertices can be found via `V(network)[V(network)$Network..Is.root]`
These are not necessarily the same as the roots of a tree/subnetwork since 
"Fake" vertices are not considered true roots for sequencing purposes as
there is no settlement associated with a "Fake" vertex.

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

The accumulate and sequence functions can be customized by providing  
`accumulator` and `selector` function definitions respectively.  These are
members of the `sequence_model` list parameter to the `sequence_plan_far` 
function.  There are several predefined functions in the R/sequence_models.R 
file.  

You can write out the NetworkPlan as nodes and segments via the `write.Network`
function.  This outputs the nodes as a csv and segments as a shapefile.  You 
can also download a scenario and all its files for analysis via the 
`download_scenario` function.   

