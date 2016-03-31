\encoding{ISO-8859-2} 
\name{Eukaryote_photo}
\alias{Eukaryote_photo}
\docType{data}
\title{An original and an augmented real genomic networks}
\description{Those networks are composed of 145 archaea, 181 eubacteria with the addition of 19 photosynthetic eukaryotes to illustrate the functions of this package.}
\usage{data(Eukaryote_photo)}
\format{
  A data frame with the following 2 variables:
   \describe{
    \item{g1}{An original \emph{igraph} network}
    \item{g2}{An augmented \emph{igraph} network}
  }
}

\examples{
\dontrun{
data(Eukaryote_photo)
info_network(g1,g2)
}
## Network characteristics:
## Total of new nodes in network Y: 19 
## Number of edges in network Y: 2014 
## Number of nodes in network Y: 345 
## Number of nodes in network X: 326 
## Total of pathways to investigate: 52975 
## Clustering coefficient network Y: 0.866758 
## Clustering coefficient network X: 0.8643119 
## Average degree + std in network Y: 11.67536 + 9.600091 
## Average degree + std in network X: 10.45399 + 8.502816 
## Average path length in network Y: 5.972801 
## Average path length in network X: 4.915441 
## Number of clusters in network Y: 21 
## Number of clusters in network X: 22 
## Average cluster size + std in network Y: 16.42857 + 32.51857 
## Average cluster size + std in network X: 14.81818 + 25.1711 
## Nodes distribution in network Y (first row taxa, second row count):
##        archaea      eubacteria eukaryote_photo 
##            145             181              19 
}
\keyword{datasets}
