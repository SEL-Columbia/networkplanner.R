\name{create_networkplan}
\alias{create_networkplan}
\title{Create a NetworkPlan from the components of a scenario}
\usage{
create_networkplan(metrics_df, proposed_network_df, proj4string)
}
\arguments{
  \item{metrics_df}{a dataframe corresponding to the
  metrics-local output of a scenario}

  \item{proposed_network_df}{SpatialLinesDataFrame
  corresponding to the proposed network output of a
  scenario}

  \item{proj4string}{string representing the projection (in
  proj4 format)}
}
\value{
A NetworkPlan object
}
\description{
Create a NetworkPlan from the components of a scenario
}

