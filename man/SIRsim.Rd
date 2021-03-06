% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_funcs.R
\name{SIRsim}
\alias{SIRsim}
\title{SIR simulation}
\usage{
SIRsim(t, startPops, params)
}
\arguments{
\item{t}{vector of time points to solve ODEs over}

\item{startPops}{vector of initial S, I and R sizes}

\item{params}{vector containing values for beta and gamma}
}
\value{
the 4 column data frame of SIR dynamics
}
\description{
Runs a deterministic SIR model with given times, starting population sizes and beta/gamma parameters
}

