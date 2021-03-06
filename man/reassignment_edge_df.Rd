\name{reassignment_edge_df}
\alias{reassignment_edge_df}
\title{Reassign fields from one undirected graph to another directed graph
via edge dataframes and a vertex mapping.}
\usage{
reassignment_edge_df(vertex_id_map, orig_edge_df, network, field)
}
\arguments{
  \item{vertex_id_map}{maps the original graph vertex ids
  to the new graph}

  \item{orig_edge_df}{the original graph edge dataframe}

  \item{field}{to be reassigned}
}
\value{
An edge dataframe with from/to fields corresponding to
edges in the graph to be updated along with the fields to
be added (can be used via the igraph index operator to
reassign fields to the igraph)
}
\description{
Reassign fields from one undirected graph to another
directed graph via edge dataframes and a vertex mapping.
}

