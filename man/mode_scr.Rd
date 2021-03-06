% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scr_mode.R
\name{mode_scr}
\alias{mode_scr}
\title{Computes the Posterior Mode}
\usage{
mode_scr(Data, prior, mixture)
}
\arguments{
\item{Data:}{data as in \code{\link{glwfm}}.}

\item{prior:}{priori distributions as in \code{\link{glwfm}}.}

\item{mixture:}{indicates the finite mixture to be fitted to the data. Options are: 'glw', 'gl', 'gw', 'lw', 'g', 'l', 'w'}
}
\value{
Returns a list containing the mode of the posterior distribution of the GLW finite mixture regression model, 
     the log-posterior maximum, the convergence status (as in optim) and a message (as in optim).
}
\description{
\code{mode_scr} computes the posterior mode of the GLW FM Standard Cure Rate model
}

