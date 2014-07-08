\name{create_networkplan_from_graph}
\alias{create_networkplan_from_graph}
\title{Make a NetworkPlan out of a graph and related vertex datafram}
\usage{
create_networkplan_from_graph(g, vertex_df, proj)
}
\arguments{
  \item{g}{an igraph}

  \item{vertex_df}{a dataframe of vertices with indices
  corresponding to the from/to attributes of the edges of
  the igraph}

  \item{proj}{the proj4 projection as a string}
}
\value{
A NetworkPlan
}
\description{
Make a NetworkPlan out of a graph and related vertex
datafram
}

