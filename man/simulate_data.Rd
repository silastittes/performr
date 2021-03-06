% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/performance_functions.R
\name{simulate_data}
\alias{simulate_data}
\title{Generates simulated data that can be used for practice and model validation}
\usage{
simulate_data(n_spp = 2, n_axis = 20, n_reps = 2, q_end = 0.9,
  mu_sig = 1, shape1_pr_mu = 2, shape2_pr_mu = 2, stretch_pr_mu = 2,
  min_pr_mu = -5, max_pr_mu = 5, nu_pr_shape = 8, nu_pr_scale = 3,
  pr_sig = 1)
}
\arguments{
\item{n_spp}{The number of groups (i.e. species) to be simulated}

\item{n_axis}{The number of sampling locations along environmental axis}

\item{n_reps}{The number of times each location along the environmental axis should be sampled}

\item{q_end}{Quantile [0,1] to determine how close to the ends of the axis sampling should start.}
}
\value{
A list containing the true parameters for each species (named true_params), and tidy data frame, containing columns for environmental axis, response trait, mean response (handy for plotting the curve), and groups (named sim_data).
}
\description{
Generates correctly formatted input data for the performance_stan() function.
}
\details{
All other arguments are hyper prior values used to generate species-level parameters. The arguments for the two shape parameters should remain above two, stretch above zero.
}
\examples{
simulate_data()
}
