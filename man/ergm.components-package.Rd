\name{ergm.components-package}
\alias{ergm.components-package}
\alias{ergm.components}
\docType{package}
\title{
ERGM Terms Based on Connectivity and Components
}
\description{
This package contains additional terms for the \pkg{ergm} package within the \pkg{statnet} suite of packages for network analysis.  It requires \pkg{ergm} and \pkg{network} to be used.

The terms in this package compute properties related to components and connectivity (e.g., number of components, non-central moments of component sizes, cutedges).  They can be used with the \link{ergm} \code{summary} function to calculate network descriptives, although for many of these properties it is (possibly much) more efficient to use functions in the \pkg{sna} library for this purpose.  A more typical use case is employing them for simulation and/or inference for exponential family random graph models (ERGMs); using these terms allows one to infer and simulate the effects of biases that impact component sizes, the prevalence of bridges, and other non-local properties.

This package was implemented using the \pkg{ergm.userterms} framework (Hunter et al., 2010).

We ask that this package be cited when using its terms in published results; this also helps readers reproduce the associated results.  \code{citation("ergm.components")} contains the relevant information.
}
\references{
Hunter DR, Goodreau SM, Handcock MS (2013).
{\pkg{ergm.userterms}: {A} Template Package for 
Extending \pkg{statnet}}, Journal of Statistical
Software 52(2), 1-25,
URL http://www.jstatsoft.org/v52/i02/.
}
\seealso{statnet, network, ergm, ergm-terms, sna, components, bridges}
%\examples{
%}
\keyword{ package }
\keyword{ models }
