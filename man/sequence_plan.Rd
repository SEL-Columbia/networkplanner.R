\name{sequence_plan}
\alias{sequence_plan}
\title{Sequence a NetworkPlan via breadth-first search from roots
and a selector function (to select from the "frontier"
of SpatialPoints)}
\usage{
sequence_plan(np, selector = default_selector)
}
\arguments{
  \item{np}{a NetworkPlan TODO: Remove roots once we set
  them in read method}

  \item{roots}{the indices of root vertices to sequence
  from}

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

