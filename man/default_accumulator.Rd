\name{default_accumulator}
\alias{default_accumulator}
\title{Default accumulator for accumulate method
This calculates the count of all downstream nodes}
\usage{
default_accumulator(node_df, edge_df, g, vid)
}
\arguments{
  \item{node_df}{downstream node dataframe}

  \item{edge_df}{downstream edge dataframe}

  \item{g}{the graph}

  \item{vid}{of the current vertex in the graph}
}
\value{
subset of the node_df
}
\description{
Default accumulator for accumulate method This calculates
the count of all downstream nodes
}

