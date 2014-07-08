\name{download_scenario}
\alias{download_scenario}
\title{Download scenario from networkplanner into the given directory}
\usage{
download_scenario(scenario_number, directory_name = NULL, username = NULL,
  password = NULL, np_url = "http://networkplanner.modilabs.org/")
}
\arguments{
  \item{scenario_number}{scenario_number on
  http://networkplanner.modilabs.org.}

  \item{directory_name}{path to write the downloaded
  scenario into. By default, and if directory_name==NULL,
  directory_name will be scenario number within current
  working directory.  for Windows user you can either use
  "\\" or "/" as the path separator, but a single
  back-slash "\" won't work in R environment as the path
  construct.  e.g. "C:\\User\\xxx\\" and "User/xxx/" are
  legit file path.}

  \item{username}{username to login to
  http://networkplanner.modilabs.org. If NULL, we assume
  that the scenario is public.}

  \item{password}{password associated with the previous
  username. If NULL, we assume that the scenario is
  public.}

  \item{np_url}{URL of the network planner instance to
  download scenario from. By default, it is
  http://networkplanner.modilabs.org.}
}
\description{
Download scenario from networkplanner into the given
directory
}

