\encoding{ISO-8859-2} 
\name{networkY}
\alias{networkY}
\docType{data}
\title{An augmented network Y with 15 nodes and 18 edges}
\description{This file contains one sample augmented network Y as described in: \cr 
			 Lord et al. (2016). BRIDES: a new fast algorithm and software for characterizing evolving similarity networks using \cr
			 breakthroughs, roadblocks, impasses, detours, equals and shortcuts. PloS one, 11(8), e0161474.}
\usage{data(networkX)}
\format{
  An igraph object containing network y
  \describe{
    \item{networkY}{An augmented \emph{igraph} network}
  }
}

\examples{
\dontrun{
	data(networkX)
	data(networkY)
	plot_network(networkX,networkY)
	BRIDES(networkX,networkY)
}
}
\keyword{dataset}
