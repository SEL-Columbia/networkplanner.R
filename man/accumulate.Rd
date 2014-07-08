\name{accumulate}
\alias{accumulate}
\title{Accumulate "downstream" nodal values into a summary value
for "this" node.}
\usage{
accumulate(np, accumulator = default_accumulator, validate = F)
}
\arguments{
  \item{np}{a NetworkPlan}

  \item{accumulator}{function that summarizes or combines
  all downstream node values into a new row dataframe
  attributed to "this" node the accumulator function takes
  node and edge dataframes along with the graph and vertex
  id of the current node}
}
\value{
A NetworkPlan whose network vertices and their associated
upstream edges have new attributes/values as defined by
accumulator
}
\description{
Accumulate "downstream" nodal values into a summary value
for "this" node.
}

