\name{sequence}
\alias{sequence}
\title{Sequence a NetworkPlan via breadth-first search from roots
and a selector function (to select from the "frontier"
of SpatialPoints)}
\usage{
sequence(np, roots, selector = default_selector)
}
\arguments{
  \item{np}{a NetworkPlan}

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

