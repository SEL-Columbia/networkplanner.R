\name{download_scenario}
\alias{download_scenario}
\title{Download scenario from networkplanner into the given directory}
\usage{
download_scenario(scenario_number, directory_name = NULL, userpwd = NULL,
  np_url = "http://networkplanner.modilabs.org/")
}
\arguments{
  \item{scenario_number}{scenario_number on
  http://networkplanner.modilabs.org}

  \item{directory_name}{path to write the downloaded
  scenario into. By default, and if directory_name==NULL,
  directory_name will be scenario number within current
  working directory}

  \item{userpwd}{username and password, separated by : (ex.
  USER:PASSWORD). If NULL, we assume that the scenario is
  public}

  \item{np_url}{URL of the network planner instance to
  download scenario from. By default, it is
  http://networkplanner.modilabs.org}
}
\description{
Download scenario from networkplanner into the given
directory
}

