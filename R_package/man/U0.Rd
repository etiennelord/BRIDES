\encoding{ISO-8859-2} 
\name{U0}
\alias{U0}
\docType{data}
\title{A sample directed and weighted network with node attributes }
\description{This dataset of a sample network with 8 nodes with attributes, and 9 directed and weighted edges.}
\usage{data(U0)}
\format{
  An igraph object containing network U0
  \describe{
    \item{U0}{An augmented \emph{igraph} network. Nodes attributes are 1 for the original network, and 2 for the augmented network.}
  }
}

\examples{
\dontrun{
	data(U0)
	plot_network(U0,attributes="2", legend_position="topleft")
	BRIDES(U0, attributes="2")
	
}
}
\keyword{dataset}
