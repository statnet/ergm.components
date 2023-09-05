\name{NewarkCommand}
\Rdversion{1.1}
\alias{NewarkCommand}
\docType{data}
\title{
Time-aggregated Communication Network from the WTC Disaster
}
\description{
A time-aggregated communication network from the collection coded by Butts et al.  (2007), containing interactions among responders to the 2001 World Trade Center disaster associated with Newark Airport.
}
\usage{data(NewarkCommand)}
\format{
 A single \code{\link[network]{network}} object (\code{NewarkCommand}) containing a directed communication network.  
}
\details{
Butts et al. (2007) coded transcripts of radio communications among 17 groups of responders to the World Trade Center Disaster of September 11, 2001, in which al-Qaeda terrorists hijacked four commercial airliners in the United States (crashing two into the World Trade Center complex in New York City, and a third into the Pentagon in Washington, D.C.).  This data set contains one time-aggregated network from the full set (see Butts et al. (2023) for the full data, and Renshaw et al. (2023) for additional details).  The network was constructed from all communications on the Newark Command channel, representing interactions among one group of responders.  Communications on this channel were used to secure Newark Airport (which was closed following the hijackings) and to mobilize police and fire equipment to be deployed at the WTC complex.  Each vertex in the network represents a single named communicant, with a tie from vertex \eqn{i} to vertex \eqn{j} if \eqn{i} called \eqn{j} at some point during the communication period.

The \code{NewarkCommand} network has two vertex attributes.  \code{has.valid} indicates whether the associated individual has valid communication events in the data set; in some cases, an individual may appear without a valid event because they were named in a communication that was partially unintelligible (i.e., both the sender and receiver could not be identified).  These (three) individuals are included by default, but users may wish to exclude these vertices in some cases.  \code{inst.coord} indicates whether the associated individual occupies an \emph{institutionalized coordinative role}, i.e. a formal organizational role such as a manager, dispatcher, or officer implying responsibility for coordinating the activities of others.

The edge attribute \code{calls} contains the number of radio calls observed between each (ordered) communicating pair during the event.  In this particular network, this ranges from 1 to 8 (with most edges involving 1-3 calling events). 

For further details on the data set, see Butts et al., (2007, 2023).
}
\source{
Butts C, Petrescu-Prahova M, Cross B (2023). \dQuote{Data on Radio
  Communication Dynamics in the World Trade Center Disaster.} R package
  version 1.0, \url{https://www.ncasd.org}.
}
\references{
Renshaw S, Livas S, Butts C, Petrescu-Prahova M (2023). \dQuote{Modeling
  Complex Interactions in a Disrupted Environment: Relational Events in
  the WTC Response.} \emph{Network Science}. \doi{10.1017/nws.2023.4}

  Butts C, Petrescu-Prahova M, Cross B (2007). \dQuote{Responder Communication
  Networks in the World Trade Center Disaster: Implications for
  Modeling of Communication Within Emergency Settings.} \emph{Journal of
  Mathematical Sociology}, 31(2), 121-147.

}
\examples{
#Load the data set
data(NewarkCommand)

#Plot it
plot(NewarkCommand, vertex.col = "inst.coord")
}
\keyword{datasets}
