\name{directed_networkplan}
\alias{directed_networkplan}
\title{Take an undirected NetworkPlan and make it "directed"
such that all components are trees with directed edges
to all vertices from the root}
\usage{
directed_networkplan(np, subnet_root_selector = default_root_selector)
}
\arguments{
  \item{np}{A NetworkPlan}
}
\value{
A new "directed" NetworkPlan object
}
\description{
Take an undirected NetworkPlan and make it "directed" such
that all components are trees with directed edges to all
vertices from the root
}

