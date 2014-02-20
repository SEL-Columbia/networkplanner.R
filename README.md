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
`nodes` (a SpatialPointsDataFrame)
`network` (an igraph object)
`pre-existing network` (a SpatialLinesDataFrame)

The `nodes` data.frame has the following columns, at least:
`id`
`is_root`

Functions
---
 * read_networkplan(metrics_local_csvfile, network_shapefile)
   * metrics_local_csvfile: is networkplanner's metrics-local.csv output for a given scenario
   * network_shapefile: is networkplanner's network.shp output for that same scenario
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

Usage Example
---
```
# Download scenario (eg. 531) from networkplanner
np <- read_networkplan('531/metrics-local.csv', '531/metrics.shp')
np_sequenced <- sequence_ratio(np, numerator='annualSales', denominator='Investment', sight='near')
# equivalent to:
np_sequenced <- sequence(np, function(x) { x$annualSales / x$Investment }, sight='near')
```
