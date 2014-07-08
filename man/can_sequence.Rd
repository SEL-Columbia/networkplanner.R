\name{can_sequence}
\alias{can_sequence}
\title{Determine whether a NetworkPlan meets sequencing criteria.
This means that there are no paths between "fake nodes" and
that there are no cycles.}
\usage{
can_sequence(np)
}
\arguments{
  \item{np}{a NetworkPlan}
}
\value{
A boolean indicating whether the NetworkPlan can be
sequenced
}
\description{
Determine whether a NetworkPlan meets sequencing criteria.
This means that there are no paths between "fake nodes" and
that there are no cycles.
}

