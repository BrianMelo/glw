% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ptcr_cpo.R
\name{cpo_ptcr}
\alias{cpo_ptcr}
\title{CPO of the GLW FM Promotion Time Cure Rate Model}
\usage{
cpo_ptcr(parm, Data, mixture)
}
\arguments{
\item{parm:}{sample of posterior parameters.}

\item{Data:}{the same data as in \code{\link{glw_scr}}.}

\item{mixture:}{indicates the finite mixture to be fitted to the data. Options are: 'glw', 'gl', 'gw', 'lw'.}
}
\value{
Returns the Conditional Predictive Ordinate of the GLW finite mixture model for every observation.
}
\description{
\code{cpo_ptcr} estimates the CPO based on a MCMC sample of the posterior distribuiton.
}

