\encoding{utf-8}
\name{export_network}
\alias{export_network}
\title{is an helper function to export to a text file a network}
\description{This function export to a file the given network which is compatible with Cytoscape or Gephi.
             One particularity is that it handle single node (disconnected) in the network.}
\usage{export_network(g, file, attributes = FALSE)}
\arguments{ .
  \item{g}{the igraph object network to export}
  \item{file}{the name of the file to save the network}
  \item{attributes}{export the igraph vertex tax attribute (i.e. V(g)$tax)), default: FALSE}
}

\examples{	
	\dontrun{
		data(networkX)
		export_network(networkX,file="mynetworkX.txt", attributes =TRUE)
	}
}

\keyword{export function}

