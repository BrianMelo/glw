% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glw_cpo.R
\name{cpo_glw}
\alias{cpo_glw}
\title{CPO of the GLW Finite Mixture Model}
\usage{
cpo_glw(parm, Data, mixture)
}
\arguments{
\item{parm:}{sample of posterior parameters.}

\item{Data:}{the same data as in \code{\link{glwfm}}.}

\item{mixture:}{indicates the finite mixture to be fitted to the data. Options are: 'glw', 'gl', 'gw', 'lw', 'g', 'l', 'w'.}
}
\value{
Returns the Conditional Predictive Ordinate of the GLW finite mixture model for every observation.
}
\description{
\code{cpo_glw} estimates the CPO based on a MCMC sample of the posterior distribuiton.
}

