\name{clean_networkplan}
\alias{clean_networkplan}
\title{Take an undirected NetworkPlan and return one that
guarantees that fake nodes belong to disjoint components
and that no cycles exist}
\usage{
clean_networkplan(np)
}
\arguments{
  \item{np}{A NetworkPlan}
}
\value{
A new "cleaned" NetworkPlan object
}
\description{
Take an undirected NetworkPlan and return one that
guarantees that fake nodes belong to disjoint components
and that no cycles exist
}

