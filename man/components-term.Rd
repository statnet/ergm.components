\name{components-term}
\alias{components}
\alias{dimers}
\alias{compsizesum}
\alias{InitErgmTerm.dimers}
\alias{InitErgmTerm.components}
\alias{InitErgmTerm.compsizesum}
\docType{methods}
\title{Component Size/Count Terms for the 'ergm' Package}
\description{
  These terms compute properties related to component counts and/or sizes.  In the directed case, they compute properties of weak components.
}
%\usage{
%terms.ergm()
\section{Component Size/Count Terms}{
 \describe{
  \item{\code{components()}}{\emph{Number of Components:}
    This term takes no arguments, and adds one statistic to the model that is equal to the number of (weak) components.  This term can be used with directed or undirected networks.
   }
  \item{\code{compsizesum(log=FALSE, power=2)}}{\emph{Sums of Functions of Component Sizes:}
    This term adds one statistic to the model that is equal to the sum
    \eqn{sum_c f(|V(c)|)}{\sum_c f(|V(c)|)},
    where the sum is over the (weak) components of the graph, and \eqn{|V(c)|} is the number of vertices in component \eqn{c}.  Currently, the function \eqn{f} is of the form \eqn{f(x)=x^k} if \code{log=FALSE}, or \eqn{f(x)=log(x)^k}{f(x)=\log(x)^k} if \code{log=TRUE}.  The exponent \eqn{k} is given by the \code{power} argument, which defaults to 2.  Important special cases with \code{log=FALSE} include \code{power=0}, which returns the number of components, and \code{power=2}, which returns the sum of squared component sizes; more generally, this function is proportional to the \eqn{k}th non-central moment of the component sizes.  (But note that \code{power=1} reduces to the number of vertices, and is unlikely to be useful.)  This term can be used with directed or undirected networks.
   }
  \item{\code{dimers()}}{\emph{Components of Order 2:}
    This term takes no arguments, and adds one statistic to the model that is equal to the number of (weak) components of order 2.  (Compare with \code{\link{isolates}}.)  This term can be used with directed or undirected networks.
   }
  }
}
%}
%\arguments{
%\item{\dots}{Specific to the model term.}
%}

%\details{
%}
%\value{\code{\link{ergm}} returns an object of class 'ergm' that is a list.
%}
\references{
}
\seealso{ergm, ergm-terms, isolates, sna::components, sna::component.dist}
\examples{
\dontrun{
data(florentine) #Load the Florentine families data

summary(flobusiness~components)   #Number of bridging edges in the 
                                  #business network
summary(flobusiness~dimers)       #Number of dimers
summary(flobusiness~compsizesum)  #Sum of squared component sizes

#Model with effects for components
fit <- ergm(flobusiness~edges+components) 
summary(fit)

#Model with effects for components and squared sizes
fit <- ergm(flobusiness~edges+components+compsizesum) 
summary(fit)

}
}
\keyword{models}
