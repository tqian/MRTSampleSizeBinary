% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_m_sigma.R
\name{compute_m_sigma}
\alias{compute_m_sigma}
\title{Computes "M" and "Sigma" matrices for the sandwich estimator of
variance-covariance matrix.}
\usage{
compute_m_sigma(avail_pattern, f_t, g_t, beta, alpha, p_t)
}
\arguments{
\item{avail_pattern}{A vector of length T that is the average availability at
each time point}

\item{f_t}{Defines marginal excursion effect MEE(t) under
alternative together with beta. Assumed to be matrix of size T*p.}

\item{g_t}{Defines success probability null curve together with
alpha. Assumed to be matrix of size T*q.}

\item{beta}{Length p vector that defines marginal excursion effect
MEE(t) under alternative together with f_t.}

\item{alpha}{Length q vector that defines success probability null
curve together with g_t.}

\item{p_t}{Length T vector of randomization probabilities at each time point}
}
\value{
List containing two matrices. The first is the M matrix
  and the second is the Sigma matrix.
}
\description{
A helper function for
mrt_binary_power() and mrt_binary_ss().
}
\examples{
           compute_m_sigma(tau_t_1, f_t_1, g_t_1, beta_1, alpha_1,
                                      p_t_1)
}
