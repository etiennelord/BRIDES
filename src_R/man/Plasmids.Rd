\encoding{ISO-8859-2} 
\name{Plasmids}
\alias{Plasmids}
\docType{data}
\title{An original and an augmented real genomic networks}
\description{These networks include 152 archaea, 217 eubacteria with the addition of 3 479 plasmids to illustrate the application of the SDDE package}
\usage{data(Plasmids)}
\format{
  A data frame with the following 2 variables:
  \describe{
    \item{g1}{An original \emph{igraph} network}
    \item{g2}{An augmented \emph{igraph} network}
  }
}

\examples{
\dontrun{
	data(Plasmids)
	info_network(g1,g2)
}
## Network characteristics:
## Total of new nodes in network Y: 3522 
## Number of edges in network Y: 187848 
## Number of nodes in network Y: 3848 
## Number of nodes in network X: 326 
## Total of pathways to investigate: 52975 
## Clustering coefficient network Y: 0.5591526 
## Clustering coefficient network X: 0.8643119 
## Average degree + std in network Y: 97.6341 + 119.6407 
## Average degree + std in network X: 10.45399 + 8.502816 
## Average path length in network Y: 2.903052 
## Average path length in network X: 4.915441 
## Number of clusters in network Y: 74 
## Number of clusters in network X: 22 
## Average cluster size + std in network Y: 52 + 417.644 
## Average cluster size + std in network X: 14.81818 + 25.1711 
## Nodes distribution in network Y (first row taxa, second row count):
##   archaea eubacteria    plasmid 
##       152        217       3479 
}
\keyword{dataset}
