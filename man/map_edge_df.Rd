\name{map_edge_df}
\alias{map_edge_df}
\title{Map edge dataframe fields to the attributes of a graph via a vertex id mapping}
\usage{
map_edge_df(from_edge_df, to_graph, vertex_id_map, multi_dir = TRUE)
}
\arguments{
  \item{from_edge_df}{the original graph edge dataframe}

  \item{to_graph}{the igraph that we're mapping to}

  \item{vertex_id_map}{maps the original graph vertex ids
  to the new graph (the index of the vector represents the
  original vertex id and the value at the index represents
  the new vertex id)}
}
\value{
An edge dataframe with from/to fields corresponding to
edges in the graph to be updated along with the fields to
be added (can be used via the igraph index operator to
reassign fields to the igraph)
}
\description{
Map edge dataframe fields to the attributes of a graph via
a vertex id mapping
}

