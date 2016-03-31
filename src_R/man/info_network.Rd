\encoding{ISO-8859-2}
\name{info_network}
\alias{info_network}
\title{returns additional information regarding the networks X and Y (original and augmented).}
\description{This function returns some additional information such as number of nodes and edges in the networks X and Y.}
\usage{info_network(g1,g2, taxnames)}
\arguments{
	\item{g1}{the original network X}
    \item{g2}{the augmented network Y with additional nodes (all the original nodes from X must be present in the augmented network Y)}
    \item{taxnames}{the taxon name of the nodes added to the original graph (default: we select all nodes that are not in g1)}
}

\value{Additional information such as: number of nodes in networks X and Y, clustering coefficient, average degree, number of clusters and node distributions.}


\examples{
   ## Report information about the networks of Sample_1   
	data(Sample_1)
	info_network(g1,g2)	
   ## Network characteristics:
   ## Total of new nodes in network Y: 3 
   ## Number of edges in network Y: 17 
   ## Number of nodes in network Y: 14 
   ## Number of nodes in network X: 11 
   ## Total of pathways to investigate: 55 
   ## Clustering coefficient network Y: 0 
   ## Clustering coefficient network X: 0 
   ## Average degree + std in network Y: 2.428571 + 1.01635 
   ## Average degree + std in network X: 2 + 0.8944272 
   ## Average path length in network Y: 2.119403 
   ## Average path length in network X: 2.162162 
   ## Number of clusters in network Y: 2 
   ## Number of clusters in network X: 2 
   ## Average cluster size + std in network Y: 7 + 7.071068 
   ## Average cluster size + std in network X: 5.5 + 4.949747 
   ## Nodes distribution in network Y (first row taxa, second row count):
   ##  1  2 
   ## 11  3 
}

%\keyword{ }

