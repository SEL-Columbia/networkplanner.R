\name{map_edge_attributes}
\alias{map_edge_attributes}
\title{assign edge attributes from one graph to another via a vertex id mapping}
\usage{
map_edge_attributes(from_graph, to_graph, id_map, multi_dir = TRUE)
}
\arguments{
  \item{from_graph}{the igraph to take edge attributes
  from}

  \item{to_graph}{the igraph to assign edge attributes to}

  \item{id_map}{the mapping of from_graph to to_graph
  vertex ids}

  \item{multi_dir}{whether to apply from_graph edge
  attributes in both directions to the to_graph}
}
\value{
a new graph matching the topology of the to_graph with
attributes of from_graph
}
\description{
assign edge attributes from one graph to another via a
vertex id mapping
}

