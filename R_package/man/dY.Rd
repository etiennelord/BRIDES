\encoding{ISO-8859-2} 
\name{dY}
\alias{dY}
\docType{data}
\title{An augmented directed network Y with 5 nodes and 7 edges}
\description{This dataset is composed of a directed augmented sample network with one additional node}
\usage{data(networkX)}
\format{
  An igraph object containing network Y
  \describe{
    \item{dY}{An augmented \emph{igraph} network}
  }
}

\examples{
\dontrun{
	data(dX)
	data(dY)
	plot_network(dX,dY, legend_position="topleft")
	BRIDES(dX,dY)
	# (see dX for results)
}
}
\keyword{dataset}
