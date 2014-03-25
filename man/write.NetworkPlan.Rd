\name{write.NetworkPlan}
\alias{write.NetworkPlan}
\title{Write line shapfile from networkplan object into the given directory}
\usage{
write.NetworkPlan(np, directory_name, nodeFormat = "csv",
  edgeFormat = "shp", includeFake = FALSE)
}
\arguments{
  \item{np}{a NetworkPlan}

  \item{directory_name}{path to write the downloaded
  scenario into. By default, and if directory_name==NULL,
  directory_name will be scenario number within current
  working directory}

  \item{nodeFormat}{a string indicating the type of output
  file for the nodes in NetworkPlan, only support 'csv' for
  now.}

  \item{edgeFormat}{a string indicating the type of output
  file for the edges in NetworkPlan, only support 'shp' for
  now.}

  \item{includeFake}{a boolean indicating whether to output
  the fake node in the vertex/node csv file, Default is set
  to False}
}
\description{
Write line shapfile from networkplan object into the given
directory
}

