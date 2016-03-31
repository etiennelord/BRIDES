\encoding{ISO-8859-2} 
\name{Viruses}
\alias{Viruses}
\docType{data}
\title{An original and an augmented real genomic networks}
\description{These networks include 145 archaea, 181 eubacteria with the addition of 1658 viruses to illustrate the application of the SDDE package.}
\usage{data(Viruses)}
\format{
  A data frame with the following 2 variables.
  \describe{
    \item{g1}{An original \emph{igraph} network}
    \item{g2}{An augmented \emph{igraph} network}
  }
}

\examples{
\dontrun{
	data(Viruses)
	info_network(g1,g2)
}
## Network characteristics:
## Total of new nodes in network Y: 1658 
## Number of edges in network Y: 12054 
## Number of nodes in network Y: 1984 
## Number of nodes in network X: 326 
## Total of pathways to investigate: 52975 
## Clustering coefficient network Y: 0.8010996 
## Clustering coefficient network X: 0.8643119 
## Average degree + std in network Y: 12.15121 + 15.62558 
## Average degree + std in network X: 10.45399 + 8.502816 
## Average path length in network Y: 5.233931 
## Average path length in network X: 4.915441 
## Number of clusters in network Y: 320 
## Number of clusters in network X: 22 
## Average cluster size + std in network Y: 6.2 + 20.8493 
## Average cluster size + std in network X: 14.81818 + 25.1711 
## Nodes distribution in network Y (first row taxa, second row count):
##   archaea eubacteria      virus 
##       145        181       1658 
}
\keyword{datasets}
