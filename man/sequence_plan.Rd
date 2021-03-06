\name{sequence_plan}
\alias{sequence_plan}
\title{Sequence a NetworkPlan via breadth-first search from roots
and a selector function (to select from the "frontier"
of SpatialPoints)}
\usage{
sequence_plan(np, selector = default_selector, validate = F)
}
\arguments{
  \item{np}{a NetworkPlan}

  \item{selector}{function that returns which vertex (by
  id) in the "frontier" of the search gets selected next
  based on SpatialPointsDataFrame}
}
\value{
A NetworkPlan whose nodes SpatialPointsDataFrame has a
sequence column and values
}
\description{
Sequence a NetworkPlan via breadth-first search from roots
and a selector function (to select from the "frontier" of
SpatialPoints)
}

