\name{sequence_plan_far}
\alias{sequence_plan_far}
\title{Perform a "far-sighted" sequencing by combining calls to
accumulate (to accumulate downstream demand for each vertex)
and sequence_plan each getting a custom function to perform
the accumulation/sequencing}
\usage{
sequence_plan_far(np, sequence_model = default_sequence_model, validate = T)
}
\arguments{
  \item{np}{a NetworkPlan}

  \item{accumulated_field}{name of field to set accumulated
  value}

  \item{sequence_model}{a list with 3 members: accumulator
  function to accumulate downstream values (see
  \code{accumulate}) accumulator_field name of attribute
  for accumulated value to go in vertices and edges
  selector function to perform sequencing (see
  \code{sequence_plan})}

  \item{validate}{whether to check if the networkplan can
  be sequenced}
}
\value{
A NetworkPlan whose network vertices have a sequence value
based on the rollout_functions.  The edges of the network
should also have the accumulated_field set to the
appropriate value.
}
\description{
Perform a "far-sighted" sequencing by combining calls to
accumulate (to accumulate downstream demand for each
vertex) and sequence_plan each getting a custom function to
perform the accumulation/sequencing
}
\seealso{
\code{\link{sequence_plan}}

\code{\link{accumulate}}
}

