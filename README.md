networkplanner.R
================

Library to pre- and post-process data for [Network Planner](http://networkplanner.modilabs.org) usage.


How to install the library
---
```
install.packages('devtools')
require(devtools)
instalL_github('networkplanner.R', 'SEL-Columbia')
```

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
?? branch identify
?? grid-length
?? units

Usage Example
---
```
# Download scenario (eg. 531) from networkplanner
np <- read_networkplan('531/metrics-local.csv', '531/metrics.shp')
np_sequenced <- sequence_ratio(np, numerator='annualSales', denominator='Investment', sight='near')
# equivalent to:
np_sequenced <- sequence(np, function(x) { x$annualSales / x$Investment }, sight='near')
```
