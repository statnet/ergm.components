\name{gDc}
\Rdversion{1.1}
\alias{gDc}
\docType{data}
\title{
Simulated Aggregation States of gammaD-Crystallin
}
\description{
This data set contains two network objects, each of which contains an aggregation graph from a snapshot from an atomistic simulation of aggregation in human \eqn{\gamma}{gamma}D-crystallin; the two objects respectively represent models of the wild-type and W42R variants.
}
\usage{data(gDc)}
\format{
 Two 375-node \code{\link[network]{network}} objects: \code{gDcAgg.WT} and \code{gDcAgg.W42R}
}
\details{
Wong et al. (2019) use multi-conformational Monte Carlo methods to simulate equilibrium draws from an atomistic model of the human eye lens protein \eqn{\gamma}{gamma}D-crystallin (gDc), a structural protein necessary to sustain the refractive index of the eye lens.  Irreversible aggregation of gDc leads to cataract, motivating the study of interaction among gDc monomers.  Frames from Wong et al.'s simulation were converted to \emph{aggregation graphs,} and were studied by Diessner et al. (2023) using network Hamiltonian models (i.e., ERGMs).

This data set contains two representative frames from the two sets of simulated draws.  One (\code{gDcAgg.WT}) represents the wild-type protein, while the other (\code{gDcAgg.W42R}) represents the cataract-associated W42R variant.  In each case, the network contains 375 vertices (representing protein monomers), with edges between proteins that were interacting (as determined by a proximity cutoff; see references for details).  Both simulations were performed under periodic boundary conditions at 310K, with an effective concentration of 200g/L.
}
\source{
Diessner, Elizabeth M.; Freites, J. Alfredo; Tobias, Douglas J.; and Butts, Carter T.  (2023).  \dQuote{Network Hamiltonian Models for Unstructured Protein Aggregates, with Application to \eqn{\gamma}{gamma}D-Crystallin.}  \emph{Journal of Physical Chemistry, B,} 127(3), 685-697.  \doi{10.1021/acs.jpcb.2c07672}
}
\references{
Wong, Eric K.; Prytkova, Vera; Freites, J. Alfredo; Butts, Carter T. and Tobias, Douglas J.  (2019).  \dQuote{Molecular Mechanism of Aggregation of the Cataract-related \eqn{\gamma}{gamma}D-Crystallin W24R Variant from Multiscale Atomistic Simulations.} \emph{Biochemistry}, 58(35), 3691-3699.  \doi{10.1021/acs.biochem.9b00208}
}
\examples{
#Load the data
data(gDc)

#Plot the two networks
op <- par(mfrow=c(1,2))
plot(gDcAgg.WT, main = "Wild Type")
plot(gDcAgg.W42R, main = "W42R", vertex.col = 4)
par(op)
}
\keyword{datasets}
