\name{create_graph}
\alias{create_graph}
\title{create graph from node dataframe (NOT sp) and segment_matrix}
\usage{
create_graph(metrics_df, segment_matrix, ids,
  proj4string = "+proj=longlat +datum=WGS84 +ellps=WGS84")
}
\arguments{
  \item{metrics_df}{dataframe representing the nodes (with
  attributes) of the graph}

  \item{segment_matrix}{NX2X2 matrix of segments X
  end_points X coordinates representing a network}

  \item{ids}{set of ids associated with each segment in
  matrix to be associated with edges}
}
\value{
an igraph object of merged metrics_df nodes and
segment_matrix segments
}
\description{
create graph from node dataframe (NOT sp) and
segment_matrix
}

