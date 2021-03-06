% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glw_mode.R
\name{mode_glw}
\alias{mode_glw}
\title{Computes the Posterior Mode}
\usage{
mode_glw(Data, prior, mixture)
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
\code{Mode_glw} computes the posterior mode of the GLW finite mixture model
}

