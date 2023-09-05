\name{bridges-term}
\alias{bridges}
\alias{InitErgmTerm.bridges}
\docType{methods}
\title{Bridging Edge (Cutedge) Term for the 'ergm' Package}
\description{
  This is a term for the \code{\link{ergm}} package that counts the number of edges whose removal increases the number of components (i.e., bridges or cutedges).  In the directed case, weak connectivity is employed.
}
\section{Bridging Edge Term}{
%\usage{
%terms.ergm()
 \describe{
  \item{\code{bridges(dyadic=TRUE)}}{\emph{Bridges:}
    In the undirected case, this term adds one statistic to the model equal to the number of bridging edges (aka cutedges).
    
    In the directed case, the term behavior depends on the \code{dyadic} argument.  When \code{dyadic==TRUE}, the term adds one statistic to the model equal to the number of semi-bridging edges (i.e., bridging edges in the underlying semigraph); thus, each dyad contribute either 0 or 1 bridge to the count.  When \code{dyadic==FALSE}, the term adds one statistic to the model equal to the number of weak cutedges.  These statistics are closely related: a weak cutedge is an edge within an asymmetric semi-bridging dyad (i.e., a semi-bridging edge with no reciprocating edge).  However, note that a graph may have multiple semi-bridges with no cutedges, the latter being a stronger condition.
   }
}
%}
}
%\arguments{
%\item{\dots}{Specific to the model term.}
%}

%\details{
%}
%\value{
%}
\references{
Berge, Claude.  (1966).  \emph{The Theory of Graphs.} New York: John Wiley and Sons.
}
\seealso{ergm, ergm-terms}
\examples{
\dontrun{
data(florentine) #Load the Florentine families data

summary(flomarriage~bridges)  #Number of bridging edges in the 
                              #marriage network

#Model with effects for bridging
fit <- ergm(flomarriage~edges+bridges) 
summary(fit)

#Identify the bridging edges, using ergm() changescores
cs<-ergmMPLE(flomarriage~bridges,output="array")$predictor[,,1]
diag(cs)<-0
cs[lower.tri(cs)]<-t(cs)[lower.tri(cs)]  #Matrix is upper-triangular
                                         #by default
plot(flomarriage,edge.col=1+(cs>0))      #Bridges are red
hist(cs)                                 #Some changes can impact
                                         #multiple bridges!
}
}
\keyword{models}
